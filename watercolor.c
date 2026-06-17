#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* layout */
#define CANVAS_W 620
#define CANVAS_H 900
#define PANEL_W  300
#define WIN_W (CANVAS_W + PANEL_W)
#define WIN_H CANVAS_H

#define MAX_LAYERS 6

/* sim constants */
#define FLOW_RATE        9.0f   /* how fast water equalizes across cells */
#define MAX_OUT_FRACTION 0.6f   /* max fraction of a cell's water that can leave per substep */
#define TEX_FLOW_BIAS    0.10f  /* how much paper grain perturbs flow direction */
#define ABSORB_RATE      0.45f  /* base drying/absorption rate per second */
#define FIX_RATE         0.06f  /* baseline pigment staining rate per second */
#define EPS_WATER        0.004f
#define EPS_SUS          0.0008f
#define SUBSTEPS         3
#define SIM_DT           (1.0f/60.0f)

static int   g_running = 1;
static Uint32 g_seed;

/* utility */

static inline float clampf(float v, float lo, float hi){ return v<lo?lo:(v>hi?hi:v); }
static inline float maxf(float a,float b){ return a>b?a:b; }
static inline float minf(float a,float b){ return a<b?a:b; }

static inline float hashf(int x, int y, unsigned int seed){
    unsigned int h = (unsigned int)x * 374761393u + (unsigned int)y * 668265263u + seed * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    h ^= (h >> 16);
    return (h & 0xFFFFFFu) / (float)0xFFFFFFu;
}

/* bilinear value noise at a given cell size, returns 0..1 */
static float valueNoise(float x, float y, float cellSize, unsigned int seed){
    float gx = x / cellSize, gy = y / cellSize;
    int x0 = (int)floorf(gx), y0 = (int)floorf(gy);
    float fx = gx - x0, fy = gy - y0;
    float v00 = hashf(x0,   y0,   seed);
    float v10 = hashf(x0+1, y0,   seed);
    float v01 = hashf(x0,   y0+1, seed);
    float v11 = hashf(x0+1, y0+1, seed);
    /* smoothstep interpolation for organic look */
    float sx = fx*fx*(3-2*fx);
    float sy = fy*fy*(3-2*fy);
    float a = v00 + (v10-v00)*sx;
    float b = v01 + (v11-v01)*sx;
    return a + (b-a)*sy;
}

/* paper */

static float *g_paperTex; /* CANVAS_W*CANVAS_H, 0..1 grain/absorbency field */
static SDL_Color g_paperColor = {247, 241, 224, 255};

static void generatePaper(unsigned int seed){
    if(!g_paperTex) g_paperTex = (float*)malloc(sizeof(float)*CANVAS_W*CANVAS_H);
    for(int y=0;y<CANVAS_H;y++){
        for(int x=0;x<CANVAS_W;x++){
            float n = 0.0f;
            n += valueNoise((float)x,(float)y, 90.0f, seed+1) * 0.35f;
            n += valueNoise((float)x,(float)y, 40.0f, seed+2) * 0.25f;
            n += valueNoise((float)x,(float)y, 16.0f, seed+3) * 0.22f;
            n += valueNoise((float)x,(float)y, 5.0f,  seed+4) * 0.18f;
            /* fine fibrous speckle, slightly anisotropic to mimic
               cold-press paper grain direction */
            float fiber = valueNoise((float)x*1.0f, (float)y*2.6f, 2.2f, seed+5);
            n = n*0.82f + fiber*0.18f;
            n = clampf(n, 0.0f, 1.0f);
            g_paperTex[y*CANVAS_W+x] = n;
        }
    }
}

/* Each layer is its own little fluid simulation. "sus*" is pigment mass
   currently suspended/mobile in the water; "dep*" is pigment that has
   already stained into the paper fibers and is permanent. Both are
   tracked in a simple CMY absorption space so multiple pigments mixed
   wet-on-wet combine the way real transparent paints do. 
*/

typedef struct {
    float *water;
    float *susC, *susM, *susY;
    float *depC, *depM, *depY;
    int   *activeFlag;   /* index into activeList, or -1 */
    int   *activeList;
    int    activeCount;
    int    activeCap;
    int    visible;
    int    dirtySinceComposite; /* set whenever this layer changed visually */
    char   name[24];
} Layer;

static Layer g_layers[MAX_LAYERS];
static int   g_layerCount = 0;
static int   g_activeLayer = 0;

static void layerInit(Layer *L, const char *name){
    int n = CANVAS_W*CANVAS_H;
    L->water = (float*)calloc(n, sizeof(float));
    L->susC  = (float*)calloc(n, sizeof(float));
    L->susM  = (float*)calloc(n, sizeof(float));
    L->susY  = (float*)calloc(n, sizeof(float));
    L->depC  = (float*)calloc(n, sizeof(float));
    L->depM  = (float*)calloc(n, sizeof(float));
    L->depY  = (float*)calloc(n, sizeof(float));
    L->activeFlag = (int*)malloc(sizeof(int)*n);
    for(int i=0;i<n;i++) L->activeFlag[i] = -1;
    L->activeCap = 4096;
    L->activeList = (int*)malloc(sizeof(int)*L->activeCap);
    L->activeCount = 0;
    L->visible = 1;
    L->dirtySinceComposite = 1;
    snprintf(L->name, sizeof(L->name), "%s", name);
}

static void layerFree(Layer *L){
    free(L->water); free(L->susC); free(L->susM); free(L->susY);
    free(L->depC); free(L->depM); free(L->depY);
    free(L->activeFlag); free(L->activeList);
    memset(L, 0, sizeof(Layer));
}

static inline void activateCell(Layer *L, int idx){
    if(L->activeFlag[idx] >= 0) return;
    if(L->activeCount >= L->activeCap){
        L->activeCap *= 2;
        L->activeList = (int*)realloc(L->activeList, sizeof(int)*L->activeCap);
    }
    L->activeList[L->activeCount] = idx;
    L->activeFlag[idx] = L->activeCount;
    L->activeCount++;
}

static inline void deactivateCell(Layer *L, int idx){
    int pos = L->activeFlag[idx];
    if(pos < 0) return;
    int last = L->activeList[L->activeCount-1];
    L->activeList[pos] = last;
    L->activeFlag[last] = pos;
    L->activeFlag[idx] = -1;
    L->activeCount--;
}

static int g_neigh8dx[8] = {-1,1,0,0,-1,-1,1,1};
static int g_neigh8dy[8] = {0,0,-1,1,-1,1,-1,1};
static float g_neighW[8]  = {1,1,1,1,0.7071f,0.7071f,0.7071f,0.7071f};

/* one fixed-size physics substep across a layer's currently active cells */
static void simulateLayerStep(Layer *L, float dt){
    if(L->activeCount == 0) return;

    /* snapshot active list, cells woken this step are processed
       next frame */
    static int *snapshot = NULL;
    static int snapCap = 0;
    if(snapCap < L->activeCount){
        snapCap = L->activeCount*2;
        snapshot = (int*)realloc(snapshot, sizeof(int)*snapCap);
    }
    int n = L->activeCount;
    memcpy(snapshot, L->activeList, sizeof(int)*n);

    float diffs[8];

    for(int i=0;i<n;i++){
        int idx = snapshot[i];
        if(L->activeFlag[idx] < 0) continue; /* deactivated earlier this loop */

        int x = idx % CANVAS_W;
        int y = idx / CANVAS_W;

        float w = L->water[idx];
        float susC = L->susC[idx], susM = L->susM[idx], susY = L->susY[idx];
        float tex = g_paperTex[idx];
        float effSelf = w - tex*TEX_FLOW_BIAS;

        /* flow to neighbors, proportional to downhill slope */
        if(w > EPS_WATER){
            float sumPos = 0.0f;
            int validN[8]; int nn=0;
            for(int k=0;k<8;k++){
                int nx = x + g_neigh8dx[k];
                int ny = y + g_neigh8dy[k];
                if(nx<0||nx>=CANVAS_W||ny<0||ny>=CANVAS_H) continue;
                int nidx = ny*CANVAS_W+nx;
                float effN = L->water[nidx] - g_paperTex[nidx]*TEX_FLOW_BIAS;
                float d = (effSelf - effN) * g_neighW[k];
                if(d > 0.0f){
                    diffs[nn] = d;
                    validN[nn] = nidx;
                    sumPos += d;
                    nn++;
                }
            }
            if(nn > 0 && sumPos > 1e-6f){
                float outflow = minf(w * MAX_OUT_FRACTION, FLOW_RATE * sumPos * dt);
                float concC = susC / w, concM = susM / w, concY = susY / w;
                for(int k=0;k<nn;k++){
                    float share = outflow * (diffs[k]/sumPos);
                    int nidx = validN[k];
                    L->water[nidx] += share;
                    L->susC[nidx]  += concC*share;
                    L->susM[nidx]  += concM*share;
                    L->susY[nidx]  += concY*share;
                    activateCell(L, nidx);
                }
                w    -= outflow;
                susC -= concC*outflow;
                susM -= concM*outflow;
                susY -= concY*outflow;
            }
        }

        /* absorption / evaporation / pigment fixation */
        float lossFrac = ABSORB_RATE * (0.35f + 0.75f*tex) * dt;
        float waterLoss = w * lossFrac;
        float wAfter = w - waterLoss;
        if(wAfter < 0) wAfter = 0;

        float fixFrac = FIX_RATE*dt + (w>1e-6f ? (waterLoss/w) : 1.0f);
        fixFrac = clampf(fixFrac, 0.0f, 1.0f);
        float fixedC = susC*fixFrac, fixedM = susM*fixFrac, fixedY = susY*fixFrac;
        susC -= fixedC; susM -= fixedM; susY -= fixedY;
        L->depC[idx] += fixedC;
        L->depM[idx] += fixedM;
        L->depY[idx] += fixedY;

        if(wAfter < EPS_WATER){
            /* fully dried this step -- dump any remaining suspension */
            L->depC[idx] += susC; L->depM[idx] += susM; L->depY[idx] += susY;
            susC = susM = susY = 0.0f;
            wAfter = 0.0f;
        }

        L->water[idx] = wAfter;
        L->susC[idx] = susC; L->susM[idx] = susM; L->susY[idx] = susY;

        if(wAfter < EPS_WATER && susC < EPS_SUS && susM < EPS_SUS && susY < EPS_SUS){
            deactivateCell(L, idx);
        }
        L->dirtySinceComposite = 1;
    }
}

/* brush */

static SDL_Color g_currentColor = {38, 70, 160, 255}; /* ultramarine-ish default */
static int   g_brushRadius = 16;
static float g_brushWetness = 0.65f; /* 0.1 (dry) .. 1.0 (very wet) */
static float g_brushLoad = 1.0f;     /* depletes while dragging a stroke */
static int   g_painting = 0;
static float g_lastPaintX=-1, g_lastPaintY=-1;

/* recompute pigment absorption whenever the active color changes */
static float g_curCp=0, g_curMp=0, g_curYp=0;
static void setCurrentColor(Uint8 r, Uint8 g, Uint8 b){
    g_currentColor.r=r; g_currentColor.g=g; g_currentColor.b=b;
    g_curCp = 1.0f - r/255.0f;
    g_curMp = 1.0f - g/255.0f;
    g_curYp = 1.0f - b/255.0f;
}

static void stampBrush(float cx, float cy){
    Layer *L = &g_layers[g_activeLayer];
    int R = g_brushRadius;
    float Cp = g_curCp, Mp = g_curMp, Yp = g_curYp;

    float waterAmt   = (0.55f + 0.85f*g_brushWetness) * g_brushLoad;
    float pigmentAmt = (0.42f + 0.40f*g_brushWetness) * g_brushLoad;

    int x0 = (int)floorf(cx-R), x1=(int)ceilf(cx+R);
    int y0 = (int)floorf(cy-R), y1=(int)ceilf(cy+R);
    x0 = (x0<0)?0:x0; y0=(y0<0)?0:y0;
    x1 = (x1>=CANVAS_W)?CANVAS_W-1:x1; y1=(y1>=CANVAS_H)?CANVAS_H-1:y1;

    float Rf = (float)R + 0.5f;
    for(int y=y0;y<=y1;y++){
        for(int x=x0;x<=x1;x++){
            float dx = x-cx, dy = y-cy;
            float d = sqrtf(dx*dx+dy*dy);
            if(d > Rf) continue;
            float t = d/Rf;
            float weight = 1.0f - t*t*(3-2*t); /* smoothstep falloff, soft round tip */
            if(weight <= 0.001f) continue;
            /* slight texture jitter so wet edge isn't a perfect circle */
            float jitter = 0.85f + 0.3f*g_paperTex[y*CANVAS_W+x];
            weight *= jitter;

            int idx = y*CANVAS_W+x;
            L->water[idx] += waterAmt * weight;
            L->susC[idx]  += Cp * pigmentAmt * weight;
            L->susM[idx]  += Mp * pigmentAmt * weight;
            L->susY[idx]  += Yp * pigmentAmt * weight;
            activateCell(L, idx);
        }
    }
    L->dirtySinceComposite = 1;
}

static void paintStroke(float x0,float y0,float x1,float y1){
    float dx=x1-x0, dy=y1-y0;
    float dist = sqrtf(dx*dx+dy*dy);
    float step = maxf(1.0f, g_brushRadius*0.35f);
    int steps = (int)(dist/step);
    if(steps < 1) steps = 1;
    for(int i=0;i<=steps;i++){
        float t = (float)i/steps;
        stampBrush(x0+dx*t, y0+dy*t);
    }
    /* brush slowly runs dry across a long stroke */
    g_brushLoad -= dist * 0.0010f;
    g_brushLoad = clampf(g_brushLoad, 0.12f, 1.0f);
}

/* compositing */

static SDL_Texture *g_canvasTex = NULL;
static Uint32 *g_pixels = NULL;
static int g_canvasDirty = 1;

static void compositeCanvas(void){
    for(int y=0;y<CANVAS_H;y++){
        for(int x=0;x<CANVAS_W;x++){
            int idx = y*CANVAS_W+x;
            float tex = g_paperTex[idx];
            float shade = 0.93f + 0.09f*tex; /* paper grain shading */
            float r = g_paperColor.r*shade;
            float g = g_paperColor.g*shade;
            float b = g_paperColor.b*shade;
            float wetGlow = 0.0f;

            for(int li=0; li<g_layerCount; li++){
                Layer *L = &g_layers[li];
                if(!L->visible) continue;
                float w = L->water[idx];

                /* granulation: deposited pigment settles more in paper
                   valleys (low tex) and less on ridges */
                float gran = 1.0f + 0.55f*(0.5f - tex);

                float Ctot = (L->depC[idx]*gran + L->susC[idx]);
                float Mtot = (L->depM[idx]*gran + L->susM[idx]);
                float Ytot = (L->depY[idx]*gran + L->susY[idx]);

                float Cs = 1.0f - expf(-Ctot*2.1f);
                float Ms = 1.0f - expf(-Mtot*2.1f);
                float Ys = 1.0f - expf(-Ytot*2.1f);

                r *= (1.0f - Cs);
                g *= (1.0f - Ms);
                b *= (1.0f - Ys);

                if(w > wetGlow) wetGlow = w;
            }

            /* wet sheen: damp areas look slightly lighter/glossier until they dry */
            float sheen = clampf(wetGlow*0.9f, 0.0f, 0.22f);
            r = r + (255.0f-r)*sheen*0.35f;
            g = g + (255.0f-g)*sheen*0.35f;
            b = b + (255.0f-b)*sheen*0.40f;

            Uint8 R=(Uint8)clampf(r,0,255), G=(Uint8)clampf(g,0,255), B=(Uint8)clampf(b,0,255);
            g_pixels[idx] = (0xFFu<<24) | (R<<16) | (G<<8) | B; /* ARGB8888 */
        }
    }
}

/* UI */

typedef struct { const char *name; Uint8 r,g,b; } Pigment;
static Pigment g_palette[16] = {
    {"Clear Water",     250,251,253},
    {"Cadmium Yellow",  255,200, 40},
    {"Cadmium Orange",  255,140, 40},
    {"Cadmium Red",     220, 40, 40},
    {"Alizarin Crimson",155, 30, 55},
    {"Permanent Rose",  210, 60,120},
    {"Ultramarine Blue", 40, 60,140},
    {"Cerulean Blue",    40,130,180},
    {"Phthalo Blue",     10, 60,110},
    {"Phthalo Green",    10,110, 90},
    {"Sap Green",        90,130, 40},
    {"Yellow Ochre",    200,150, 60},
    {"Burnt Sienna",    140, 75, 40},
    {"Burnt Umber",      80, 55, 40},
    {"Payne's Grey",     60, 70, 85},
    {"Lamp Black",       35, 30, 30},
};

static TTF_Font *g_fontBody = NULL;
static TTF_Font *g_fontTitle = NULL;
static TTF_Font *g_fontSmall = NULL;

static int loadFonts(void){
    const char *candidates[] = {
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
        "/Library/Fonts/Arial.ttf",
        "C:\\Windows\\Fonts\\arial.ttf",
        NULL
    };
    const char *bodyPath = NULL, *boldPath = NULL;
    for(int i=0;candidates[i];i++){
        FILE *f = fopen(candidates[i], "rb");
        if(f){
            fclose(f);
            if(!boldPath && strstr(candidates[i],"Bold")) boldPath = candidates[i];
            if(!bodyPath) bodyPath = candidates[i];
            if(boldPath && bodyPath) break;
        }
    }
    if(!boldPath) boldPath = bodyPath;
    if(!bodyPath) return 0;
    g_fontTitle = TTF_OpenFont(boldPath, 21);
    g_fontBody  = TTF_OpenFont(boldPath, 14);
    g_fontSmall = TTF_OpenFont(bodyPath, 12);
    return (g_fontBody && g_fontTitle && g_fontSmall);
}

static SDL_Renderer *g_ren = NULL;

static void drawText(TTF_Font *font, int x, int y, const char *text, SDL_Color color){
    if(!font || !text || !text[0]) return;
    SDL_Surface *surf = TTF_RenderText_Blended(font, text, color);
    if(!surf) return;
    SDL_Texture *tex = SDL_CreateTextureFromSurface(g_ren, surf);
    SDL_Rect dst = {x, y, surf->w, surf->h};
    SDL_RenderCopy(g_ren, tex, NULL, &dst);
    SDL_DestroyTexture(tex);
    SDL_FreeSurface(surf);
}
static int textWidth(TTF_Font *font, const char *text){
    if(!font) return 0;
    int w=0,h=0; TTF_SizeText(font, text, &w, &h); return w;
}

static void fillRect(int x,int y,int w,int h, Uint8 r,Uint8 g,Uint8 b,Uint8 a){
    SDL_SetRenderDrawColor(g_ren, r,g,b,a);
    SDL_Rect rc = {x,y,w,h};
    SDL_RenderFillRect(g_ren, &rc);
}
static void strokeRect(int x,int y,int w,int h, Uint8 r,Uint8 g,Uint8 b,Uint8 a){
    SDL_SetRenderDrawColor(g_ren, r,g,b,a);
    SDL_Rect rc = {x,y,w,h};
    SDL_RenderDrawRect(g_ren, &rc);
}

/* layout constants */
#define PAD 16
#define PX (CANVAS_W+PAD)
#define PWIDTH (PANEL_W-2*PAD)

#define SW_COLS 4
#define SW_W ((PWIDTH-(SW_COLS-1)*6)/SW_COLS)
#define SW_H 30
#define SW_GAP 6
#define SW_TOP 66

#define SLIDER_H 10

typedef struct { int x,y,w,h; } Rect;
static inline int inRect(Rect r, int mx, int my){
    return mx>=r.x && mx<r.x+r.w && my>=r.y && my<r.y+r.h;
}

static Rect paletteRect(int i){
    int col = i % SW_COLS, row = i / SW_COLS;
    Rect r = { PX + col*(SW_W+SW_GAP), SW_TOP + row*(SW_H+SW_GAP), SW_W, SW_H };
    return r;
}
#define PALETTE_BOTTOM (SW_TOP + 4*(SW_H+SW_GAP) - SW_GAP)

#define SLIDER_R_Y (PALETTE_BOTTOM + 46)
#define SLIDER_G_Y (SLIDER_R_Y + 26)
#define SLIDER_B_Y (SLIDER_R_Y + 52)
#define PREVIEW_Y  (SLIDER_R_Y + 80)

#define BRUSH_SEC_Y (PREVIEW_Y + 56)
#define SLIDER_SIZE_Y (BRUSH_SEC_Y + 26)
#define SLIDER_WET_Y  (BRUSH_SEC_Y + 70)

#define LAYER_SEC_Y (SLIDER_WET_Y + 46)
#define LAYER_ROW_H 28
#define LAYER_ROW_GAP 4
#define ADDBTN_Y (LAYER_SEC_Y + 20)

#define LAYER_ROWS_TOP (ADDBTN_Y + 30)

static Rect layerRowRect(int i){
    Rect r = { PX, LAYER_ROWS_TOP + i*(LAYER_ROW_H+LAYER_ROW_GAP), PWIDTH, LAYER_ROW_H };
    return r;
}

#define ACTION_BTN_H 32
#define ACTION_GAP 8
#define ACTIONS_BOTTOM (WIN_H-14)
#define ACTIONS_TOP (ACTIONS_BOTTOM - 3*ACTION_BTN_H - 2*ACTION_GAP)

static Rect actionRect(int i){
    Rect r = { PX, ACTIONS_TOP + i*(ACTION_BTN_H+ACTION_GAP), PWIDTH, ACTION_BTN_H };
    return r;
}

static Rect addLayerBtnRect(void){
    Rect r = { PX, ADDBTN_Y, PWIDTH/2 - 4, 24 };
    return r;
}
static Rect delLayerBtnRect(void){
    Rect r = { PX + PWIDTH/2 + 4, ADDBTN_Y, PWIDTH/2 - 4, 24 };
    return r;
}
static Rect sliderRect(int y){
    Rect r = { PX, y, PWIDTH, SLIDER_H };
    return r;
}

/* slider drag state */
enum { DRAG_NONE, DRAG_R, DRAG_G, DRAG_B, DRAG_SIZE, DRAG_WET };
static int g_dragMode = DRAG_NONE;

static float sliderValueFromMouse(Rect r, int mx, float lo, float hi){
    float t = (float)(mx - r.x) / (float)r.w;
    t = clampf(t, 0.0f, 1.0f);
    return lo + t*(hi-lo);
}

static void drawSlider(Rect r, float val, float lo, float hi, SDL_Color fillColor){
    fillRect(r.x, r.y, r.w, r.h, 210,202,184,255);
    float t = (val-lo)/(hi-lo);
    int fw = (int)(r.w*clampf(t,0,1));
    fillRect(r.x, r.y, fw, r.h, fillColor.r, fillColor.g, fillColor.b, 255);
    strokeRect(r.x, r.y, r.w, r.h, 120,110,95,255);
    /* handle */
    int hx = r.x + fw;
    fillRect(hx-3, r.y-3, 6, r.h+6, 70,60,50,255);
}

/* ----- flash message ----- */
static char g_flashMsg[320] = "";
static float g_flashTimer = 0.0f;
static void flash(const char *msg){
    snprintf(g_flashMsg, sizeof(g_flashMsg), "%s", msg);
    g_flashTimer = 1.8f;
}

/* layer management */
static void addLayer(void){
    if(g_layerCount >= MAX_LAYERS){ flash("Maximum layers reached"); return; }
    char nm[24];
    snprintf(nm, sizeof(nm), "Layer %d", g_layerCount+1);
    layerInit(&g_layers[g_layerCount], nm);
    g_activeLayer = g_layerCount;
    g_layerCount++;
    g_canvasDirty = 1;
    flash("Layer added");
}
static void deleteActiveLayer(void){
    if(g_layerCount <= 1){ flash("At least one layer is required"); return; }
    layerFree(&g_layers[g_activeLayer]);
    for(int i=g_activeLayer; i<g_layerCount-1; i++) g_layers[i] = g_layers[i+1];
    g_layerCount--;
    if(g_activeLayer >= g_layerCount) g_activeLayer = g_layerCount-1;
    for(int i=0;i<g_layerCount;i++) g_layers[i].dirtySinceComposite = 1;
    g_canvasDirty = 1;
    flash("Layer deleted");
}
static void clearActiveLayer(void){
    Layer *L = &g_layers[g_activeLayer];
    int n = CANVAS_W*CANVAS_H;
    memset(L->water,0,sizeof(float)*n);
    memset(L->susC,0,sizeof(float)*n); memset(L->susM,0,sizeof(float)*n); memset(L->susY,0,sizeof(float)*n);
    memset(L->depC,0,sizeof(float)*n); memset(L->depM,0,sizeof(float)*n); memset(L->depY,0,sizeof(float)*n);
    for(int i=0;i<n;i++) L->activeFlag[i] = -1;
    L->activeCount = 0;
    L->dirtySinceComposite = 1;
    g_canvasDirty = 1;
    flash("Layer cleared");
}
static void newPaper(void){
    g_seed = (unsigned int)time(NULL) ^ (unsigned int)rand();
    generatePaper(g_seed);
    for(int i=0;i<g_layerCount;i++) layerFree(&g_layers[i]);
    layerInit(&g_layers[0], "Layer 1");
    g_layerCount = 1;
    g_activeLayer = 0;
    g_canvasDirty = 1;
    flash("Fresh sheet of paper");
}
static void savePainting(void){
    SDL_Surface *surf = SDL_CreateRGBSurfaceFrom(g_pixels, CANVAS_W, CANVAS_H, 32,
        CANVAS_W*sizeof(Uint32), 0x00FF0000,0x0000FF00,0x000000FF,0xFF000000);
    char path[256];
    time_t t = time(NULL);
    snprintf(path, sizeof(path), "watercolor_painting_%ld.bmp", (long)t);
    if(surf){
        SDL_SaveBMP(surf, path);
        SDL_FreeSurface(surf);
        char msg[320];
        snprintf(msg, sizeof(msg), "Saved %s", path);
        flash(msg);
        printf("Saved painting to %s\n", path);
    }
}

/* panel rendering */
static void renderPanel(void){
    fillRect(CANVAS_W, 0, PANEL_W, WIN_H, 233,227,213,255);
    strokeRect(CANVAS_W, 0, 1, WIN_H, 170,160,145,255);

    SDL_Color dark = {58,48,40,255};
    SDL_Color mid  = {95,85,75,255};

    drawText(g_fontTitle, PX, 16, "Watercolor", dark);

    drawText(g_fontBody, PX, 44, "Pigments", mid);
    for(int i=0;i<16;i++){
        Rect r = paletteRect(i);
        Pigment *p = &g_palette[i];
        fillRect(r.x,r.y,r.w,r.h, p->r,p->g,p->b,255);
        int sel = (g_currentColor.r==p->r && g_currentColor.g==p->g && g_currentColor.b==p->b);
        if(sel) strokeRect(r.x-2,r.y-2,r.w+4,r.h+4, 40,30,20,255);
        strokeRect(r.x,r.y,r.w,r.h, 120,110,95,255);
    }

    drawText(g_fontBody, PX, PALETTE_BOTTOM+18, "Custom Mix", mid);
    SDL_Color rCol={200,60,60,255}, gCol={60,160,80,255}, bCol={70,90,210,255};
    drawSlider(sliderRect(SLIDER_R_Y), g_currentColor.r, 0,255, rCol);
    drawSlider(sliderRect(SLIDER_G_Y), g_currentColor.g, 0,255, gCol);
    drawSlider(sliderRect(SLIDER_B_Y), g_currentColor.b, 0,255, bCol);

    fillRect(PX, PREVIEW_Y, PWIDTH, 28, g_currentColor.r,g_currentColor.g,g_currentColor.b,255);
    strokeRect(PX, PREVIEW_Y, PWIDTH, 28, 120,110,95,255);

    drawText(g_fontBody, PX, BRUSH_SEC_Y, "Brush", mid);
    char buf[64];
    snprintf(buf,sizeof(buf),"Size: %d px", g_brushRadius);
    drawText(g_fontSmall, PX, SLIDER_SIZE_Y-16, buf, mid);
    drawSlider(sliderRect(SLIDER_SIZE_Y), (float)g_brushRadius, 2,70, (SDL_Color){120,110,95,255});

    snprintf(buf,sizeof(buf),"Wetness: %d%%", (int)(g_brushWetness*100));
    drawText(g_fontSmall, PX, SLIDER_WET_Y-16, buf, mid);
    drawSlider(sliderRect(SLIDER_WET_Y), g_brushWetness, 0.1f,1.0f, (SDL_Color){90,140,180,255});

    drawText(g_fontBody, PX, LAYER_SEC_Y, "Layers", mid);
    Rect addR = addLayerBtnRect(), delR = delLayerBtnRect();
    fillRect(addR.x,addR.y,addR.w,addR.h, 200,215,190,255);
    strokeRect(addR.x,addR.y,addR.w,addR.h, 120,110,95,255);
    drawText(g_fontSmall, addR.x+10, addR.y+4, "+ Add", dark);
    fillRect(delR.x,delR.y,delR.w,delR.h, 222,196,196,255);
    strokeRect(delR.x,delR.y,delR.w,delR.h, 120,110,95,255);
    drawText(g_fontSmall, delR.x+10, delR.y+4, "Delete", dark);

    for(int i=0;i<g_layerCount;i++){
        Rect r = layerRowRect(i);
        int active = (i==g_activeLayer);
        fillRect(r.x,r.y,r.w,r.h, active?(Uint8)210:(Uint8)238, active?(Uint8)225:(Uint8)233, active?(Uint8)205:(Uint8)222,255);
        strokeRect(r.x,r.y,r.w,r.h, 120,110,95,255);
        /* visibility toggle box */
        Rect eye = {r.x+4, r.y+4, 20,20};
        fillRect(eye.x,eye.y,eye.w,eye.h, g_layers[i].visible?230:170, g_layers[i].visible?230:170, g_layers[i].visible?230:170,255);
        strokeRect(eye.x,eye.y,eye.w,eye.h, 100,90,80,255);
        if(g_layers[i].visible) drawText(g_fontSmall, eye.x+5, eye.y+1, "o", dark);
        drawText(g_fontSmall, r.x+32, r.y+6, g_layers[i].name, dark);
    }

    const char *actionLabels[3] = {"New Paper", "Clear Layer", "Save Painting"};
    SDL_Color actionColors[3] = {{215,210,225,255},{225,215,195,255},{200,220,205,255}};
    for(int i=0;i<3;i++){
        Rect r = actionRect(i);
        fillRect(r.x,r.y,r.w,r.h, actionColors[i].r,actionColors[i].g,actionColors[i].b,255);
        strokeRect(r.x,r.y,r.w,r.h, 120,110,95,255);
        int tw = textWidth(g_fontBody, actionLabels[i]);
        drawText(g_fontBody, r.x + (r.w-tw)/2, r.y+7, actionLabels[i], dark);
    }

    if(g_flashTimer > 0.0f){
        SDL_Color flashCol = {70,60,50,(Uint8)clampf(g_flashTimer*180,0,255)};
        drawText(g_fontSmall, PX, ACTIONS_TOP-20, g_flashMsg, flashCol);
    }
}

/* returns 1 if the click/drag was consumed by the UI panel */
static int handlePanelMouseDown(int mx, int my){
    if(mx < CANVAS_W) return 0;

    for(int i=0;i<16;i++){
        Rect r = paletteRect(i);
        if(inRect(r,mx,my)){
            setCurrentColor(g_palette[i].r, g_palette[i].g, g_palette[i].b);
            g_brushLoad = 1.0f;
            return 1;
        }
    }
    Rect rr = sliderRect(SLIDER_R_Y), rg = sliderRect(SLIDER_G_Y), rb = sliderRect(SLIDER_B_Y);
    Rect rsz = sliderRect(SLIDER_SIZE_Y), rwet = sliderRect(SLIDER_WET_Y);
    Rect padR = {rr.x, rr.y-6, rr.w, rr.h+12};
    Rect padG = {rg.x, rg.y-6, rg.w, rg.h+12};
    Rect padB = {rb.x, rb.y-6, rb.w, rb.h+12};
    Rect padSz = {rsz.x, rsz.y-6, rsz.w, rsz.h+12};
    Rect padWet = {rwet.x, rwet.y-6, rwet.w, rwet.h+12};

    if(inRect(padR,mx,my)){ g_dragMode=DRAG_R; setCurrentColor((Uint8)sliderValueFromMouse(rr,mx,0,255), g_currentColor.g, g_currentColor.b); return 1; }
    if(inRect(padG,mx,my)){ g_dragMode=DRAG_G; setCurrentColor(g_currentColor.r, (Uint8)sliderValueFromMouse(rg,mx,0,255), g_currentColor.b); return 1; }
    if(inRect(padB,mx,my)){ g_dragMode=DRAG_B; setCurrentColor(g_currentColor.r, g_currentColor.g, (Uint8)sliderValueFromMouse(rb,mx,0,255)); return 1; }
    if(inRect(padSz,mx,my)){ g_dragMode=DRAG_SIZE; g_brushRadius=(int)sliderValueFromMouse(rsz,mx,2,70); return 1; }
    if(inRect(padWet,mx,my)){ g_dragMode=DRAG_WET; g_brushWetness=sliderValueFromMouse(rwet,mx,0.1f,1.0f); return 1; }

    if(inRect(addLayerBtnRect(),mx,my)){ addLayer(); return 1; }
    if(inRect(delLayerBtnRect(),mx,my)){ deleteActiveLayer(); return 1; }

    for(int i=0;i<g_layerCount;i++){
        Rect r = layerRowRect(i);
        if(inRect(r,mx,my)){
            Rect eye = {r.x+4, r.y+4, 20,20};
            if(inRect(eye,mx,my)){
                g_layers[i].visible = !g_layers[i].visible;
                g_layers[i].dirtySinceComposite = 1;
                g_canvasDirty = 1;
            } else {
                g_activeLayer = i;
            }
            return 1;
        }
    }
    for(int i=0;i<3;i++){
        Rect r = actionRect(i);
        if(inRect(r,mx,my)){
            if(i==0) newPaper();
            else if(i==1) clearActiveLayer();
            else savePainting();
            return 1;
        }
    }
    return 1; /* clicks anywhere in panel area are absorbed, never paint */
}

static void handlePanelMouseDrag(int mx,int my){
    (void)my;
    Rect rr = sliderRect(SLIDER_R_Y), rg = sliderRect(SLIDER_G_Y), rb = sliderRect(SLIDER_B_Y);
    Rect rsz = sliderRect(SLIDER_SIZE_Y), rwet = sliderRect(SLIDER_WET_Y);
    switch(g_dragMode){
        case DRAG_R: setCurrentColor((Uint8)sliderValueFromMouse(rr,mx,0,255), g_currentColor.g, g_currentColor.b); break;
        case DRAG_G: setCurrentColor(g_currentColor.r, (Uint8)sliderValueFromMouse(rg,mx,0,255), g_currentColor.b); break;
        case DRAG_B: setCurrentColor(g_currentColor.r, g_currentColor.g, (Uint8)sliderValueFromMouse(rb,mx,0,255)); break;
        case DRAG_SIZE: g_brushRadius=(int)sliderValueFromMouse(rsz,mx,2,70); break;
        case DRAG_WET: g_brushWetness=sliderValueFromMouse(rwet,mx,0.1f,1.0f); break;
        default: break;
    }
}

/* main */

int main(int argc, char **argv){
    (void)argc; (void)argv;
    g_seed = (unsigned int)time(NULL);
    srand(g_seed);

    if(SDL_Init(SDL_INIT_VIDEO) != 0){
        fprintf(stderr, "SDL_Init failed: %s\n", SDL_GetError());
        return 1;
    }
    if(TTF_Init() != 0){
        fprintf(stderr, "TTF_Init failed: %s\n", TTF_GetError());
    }
    if(!loadFonts()){
        fprintf(stderr, "Warning: could not find a system font; UI labels will be blank.\n"
                         "Install one (e.g. 'apt-get install fonts-dejavu-core') for full UI text.\n");
    }

    SDL_Window *win = SDL_CreateWindow("Watercolor",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        WIN_W, WIN_H, SDL_WINDOW_SHOWN);
    SDL_Renderer *ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);
    g_ren = ren;

    generatePaper(g_seed);
    layerInit(&g_layers[0], "Layer 1");
    g_layerCount = 1;
    setCurrentColor(g_currentColor.r, g_currentColor.g, g_currentColor.b);

    g_canvasTex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_ARGB8888,
        SDL_TEXTUREACCESS_STREAMING, CANVAS_W, CANVAS_H);
    g_pixels = (Uint32*)malloc(sizeof(Uint32)*CANVAS_W*CANVAS_H);

    compositeCanvas();

    Uint32 lastTicks = SDL_GetTicks();
    int frame = 0;

    while(g_running){
        SDL_Event e;
        while(SDL_PollEvent(&e)){
            if(e.type == SDL_QUIT) g_running = 0;
            else if(e.type == SDL_MOUSEBUTTONDOWN){
                if(e.button.x >= CANVAS_W){
                    handlePanelMouseDown(e.button.x, e.button.y);
                } else {
                    g_painting = 1;
                    g_brushLoad = 1.0f;
                    g_lastPaintX = e.button.x; g_lastPaintY = e.button.y;
                    stampBrush((float)e.button.x, (float)e.button.y);
                    g_canvasDirty = 1;
                }
            }
            else if(e.type == SDL_MOUSEBUTTONUP){
                g_painting = 0;
                g_dragMode = DRAG_NONE;
            }
            else if(e.type == SDL_MOUSEMOTION){
                if(g_painting){
                    paintStroke(g_lastPaintX, g_lastPaintY, (float)e.motion.x, (float)e.motion.y);
                    g_lastPaintX = e.motion.x; g_lastPaintY = e.motion.y;
                    g_canvasDirty = 1;
                } else if(g_dragMode != DRAG_NONE){
                    handlePanelMouseDrag(e.motion.x, e.motion.y);
                }
            }
            else if(e.type == SDL_MOUSEWHEEL){
                g_brushRadius += e.wheel.y*2;
                g_brushRadius = (int)clampf(g_brushRadius, 2, 70);
            }
            else if(e.type == SDL_KEYDOWN){
                SDL_Keycode k = e.key.keysym.sym;
                if(k == SDLK_ESCAPE) g_running = 0;
                else if(k == SDLK_LEFTBRACKET)  g_brushRadius = (int)clampf(g_brushRadius-2,2,70);
                else if(k == SDLK_RIGHTBRACKET) g_brushRadius = (int)clampf(g_brushRadius+2,2,70);
                else if(k == SDLK_n) newPaper();
                else if(k == SDLK_s) savePainting();
                else if(k == SDLK_c) clearActiveLayer();
                else if(k == SDLK_v){ g_layers[g_activeLayer].visible = !g_layers[g_activeLayer].visible;
                                       g_layers[g_activeLayer].dirtySinceComposite = 1; g_canvasDirty = 1; }
                else if(k >= SDLK_1 && k <= SDLK_6){
                    int idx = k - SDLK_1;
                    if(idx < g_layerCount) g_activeLayer = idx;
                }
            }
        }

        for(int li=0; li<g_layerCount; li++)
            for(int s=0;s<SUBSTEPS;s++)
                simulateLayerStep(&g_layers[li], SIM_DT/SUBSTEPS);

        for(int li=0; li<g_layerCount; li++){
            if(g_layers[li].activeCount>0 || g_layers[li].dirtySinceComposite){
                g_canvasDirty = 1;
                g_layers[li].dirtySinceComposite = 0;
            }
        }

        if(g_canvasDirty){
            compositeCanvas();
            SDL_UpdateTexture(g_canvasTex, NULL, g_pixels, CANVAS_W*sizeof(Uint32));
            g_canvasDirty = 0;
        }

        if(g_flashTimer > 0.0f) g_flashTimer -= SIM_DT;

        SDL_SetRenderDrawColor(ren, 235,230,220,255);
        SDL_RenderClear(ren);
        SDL_Rect dst = {0,0,CANVAS_W,CANVAS_H};
        SDL_RenderCopy(ren, g_canvasTex, NULL, &dst);
        renderPanel();
        SDL_RenderPresent(ren);

        frame++;
        Uint32 now = SDL_GetTicks();
        Uint32 elapsed = now-lastTicks;
        if(elapsed < 16) SDL_Delay(16-elapsed);
        lastTicks = SDL_GetTicks();

        if(frame > 5 && getenv("WC_TEST_EXIT")) g_running = 0;
    }

    SDL_DestroyTexture(g_canvasTex);
    free(g_pixels);
    for(int i=0;i<g_layerCount;i++) layerFree(&g_layers[i]);
    free(g_paperTex);
    if(g_fontBody) TTF_CloseFont(g_fontBody);
    if(g_fontTitle) TTF_CloseFont(g_fontTitle);
    if(g_fontSmall) TTF_CloseFont(g_fontSmall);
    TTF_Quit();
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}
