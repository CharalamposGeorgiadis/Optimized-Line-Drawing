#include "precomp.h"
#include "game.h"

#define LINES		1024
#define LINES_SHIFT 10
#define LINEFILE	"lines1024.dat"
#define ITERATIONS	16

#define IRand(x) ((int)(RandomFloat()*(x)))

int lx1[LINES], ly1[LINES], lx2[LINES], ly2[LINES];			// lines: start and end coordinates
int x1_, y1_, x2_, y2_;										// room for storing line backup
__int64 fitness = 0xfffffffff;								// similarity to reference image
int lidx = 0;												// current line to be mutated
float peak = 0;												// peak line rendering performance
Surface* reference, * backup;								// surfaces
BYTE* ref8;													// grayscale image for evaluation
Timer tm;													// stopwatch

BYTE* col8, * back8;

int height_1 = SCRHEIGHT - 1;
int height_8 = SCRHEIGHT - 8;
int height_16 = SCRHEIGHT - 16;
int height_24 = SCRHEIGHT - 24;
int height_32 = SCRHEIGHT - 32;
int height_33 = SCRHEIGHT - 33;
int width_1 = SCRWIDTH - 1;
int iterations_1000 = ITERATIONS * 1000;
int dims = SCRWIDTH * SCRHEIGHT;
int dims_mult = dims * 4;
int dims_div = dims / 4;


// -----------------------------------------------------------
// Mutate
// Randomly modify or replace one line.
// -----------------------------------------------------------
void MutateLine(int i)
{
    // backup the line before modifying it
    x1_ = lx1[i], y1_ = ly1[i];
    x2_ = lx2[i], y2_ = ly2[i];
    do
    {
        if (rand() & 1)
        {
            // small mutation (50% probability)
            lx1[i] += IRand(6) - 3, ly1[i] += IRand(6) - 3;
            lx2[i] += IRand(6) - 3, ly2[i] += IRand(6) - 3;
            // ensure the line stays on the screen
            lx1[i] = min(width_1, max(0, lx1[i]));
            lx2[i] = min(width_1, max(0, lx2[i]));
            ly1[i] = min(height_1, max(0, ly1[i]));
            ly2[i] = min(height_1, max(0, ly2[i]));
        }
        else
        {
            // new line (50% probability)
            lx1[i] = IRand(SCRWIDTH), lx2[i] = IRand(SCRWIDTH);
            ly1[i] = IRand(SCRHEIGHT), ly2[i] = IRand(SCRHEIGHT);
        }
    } while ((abs(lx1[i] - lx2[i]) < 3) || (abs(ly1[i] - ly2[i]) < 3));
}

void UndoMutation(int i)
{
    // restore line i to the backuped state
    lx1[i] = x1_, ly1[i] = y1_;
    lx2[i] = x2_, ly2[i] = y2_;
}

// -----------------------------------------------------------
// DrawWuLine
// Anti-aliased line rendering.
// Straight from: 
// https://www.codeproject.com/Articles/13360/Antialiasing-Wu-Algorithm
// -----------------------------------------------------------
void DrawWuLine(uint8_t* col8, int X0, int Y0, int X1, int Y1, uint8_t clrLine)
{
    /* Make sure the line runs top to bottom */
    if (Y0 > Y1)
    {
        int Temp = Y0; Y0 = Y1; Y1 = Temp;
        Temp = X0; X0 = X1; X1 = Temp;
    }

    /* Draw the initial pixel, which is always exactly intersected by
    the line and so needs no weighting */
    col8[X0 + Y0 * SCRWIDTH] = clrLine;

    /* Special-case horizontal, vertical, and diagonal lines, which
       require no weighting because they go right through the center of
       every pixel */
    int DeltaX = X1 - X0;
    int XDir = 1;

    if (DeltaX < 0)
    {
        XDir = -1;
        DeltaX = -DeltaX;
    }

    int DeltaY = Y1 - Y0;

    unsigned short ErrorAdj;
    unsigned short ErrorAccTemp, Weighting;
    /* Line is not horizontal, diagonal, or vertical */
    unsigned short ErrorAcc = 0;

    int rl = clrLine;
    int rl_shift = rl << 8;

    /* Is this an X-major or Y-major line? */
    if (DeltaY > DeltaX)
    {
        /* Y-major line; calculate 16-bit fixed-point fractional part of a
           pixel that X advances each time Y advances 1 pixel, truncating the
           result so that we won't overrun the endpoint along the X axis */
        ErrorAdj = (DeltaX << 16) / DeltaY;
        /* Draw all pixels other than the first and last */
        while (--DeltaY)
        {
            ErrorAccTemp = ErrorAcc; /* remember currrent accumulated error */
            ErrorAcc += ErrorAdj;/* calculate error for next pixel */
            if (ErrorAcc <= ErrorAccTemp)
            {
                /* The error accumulator turned over, so advance the X coord */
                X0 += XDir;
            }
            Y0++;
            /* The IntensityBits most significant bits of ErrorAcc give us the
                intensity weighting for this pixel, and the complement of the
                weighting for the paired pixel */

            Weighting = ErrorAcc >> 8;
            int offset = X0 + Y0 * SCRWIDTH;
            int offset_next = offset + XDir;

            uint8_t clrBackGround = col8[offset];
            int rr = (Weighting * (clrBackGround - rl) + rl_shift) >> 8;
            col8[offset] = rr;

            clrBackGround = col8[offset_next];
            rr = ((Weighting ^ 255) * (clrBackGround - rl) + rl_shift) >> 8;
            col8[offset_next] = rr;
        }
        /* Draw the final pixel, which is always exactly intersected by the line
        and so needs no weighting */
        col8[X1 + Y1 * SCRWIDTH] = clrLine;
        return;
    }
    /* It's an X-major line; calculate 16-bit fixed-point fractional part of a
        pixel that Y advances each time X advances 1 pixel, truncating the
        result to avoid overrunning the endpoint along the X axis */
    ErrorAdj = (DeltaY << 16) / DeltaX;
    /* Draw all pixels other than the first and last */
    while (--DeltaX)
    {
        ErrorAccTemp = ErrorAcc; /* remember currrent accumulated error */
        ErrorAcc += ErrorAdj; /* calculate error for next pixel */
        if (ErrorAcc <= ErrorAccTemp)
        {
            /* The error accumulator turned over, so advance the Y coord */
            Y0++;
        }
        /* X-major, so always advance X */
        /* The IntensityBits most significant bits of ErrorAcc give us the
        intensity weighting for this pixel, and the complement of the
        weighting for the paired pixel */
        X0 += XDir;

        Weighting = ErrorAcc >> 8;
        int offset = X0 + Y0 * SCRWIDTH;
        int offset_next = offset + SCRWIDTH;

        uint8_t clrBackGround = col8[offset];
        int rr = (Weighting * (clrBackGround - rl) + rl_shift) >> 8;
        col8[offset] = rr;

        clrBackGround = col8[offset_next];
        rr = ((Weighting ^ 255) * (clrBackGround - rl) + rl_shift) >> 8;
        col8[offset_next] = rr;
    }
    /* Draw the final pixel, which is always exactly intersected by the line
       and so needs no weighting */
    col8[X1 + Y1 * SCRWIDTH] = clrLine;
}

// -----------------------------------------------------------
// Fitness evaluation
// Compare current generation against reference image.
// -----------------------------------------------------------
__int64 Game::Evaluate()
{
    // compare to reference using SIMD magic. don't worry about it, it's fast.
    const int quads = dims_div;
    const __m128i mask4 = _mm_set1_epi32(0x000000FF);
    __m128i* A4 = (__m128i*)col8;
    __m128i* B4 = (__m128i*)ref8;
    union { __m128i diff4; int diff[4]; };
    diff4 = _mm_set1_epi8(0);

    for (int i = 0; i < quads; i++)
    {
        __m128i A8 = _mm_loadu_si128((__m128i*)(&col8[i * 4]));
        __m128i B8 = _mm_loadu_si128((__m128i*)(&ref8[i * 4]));
        __m128i d2 = _mm_abs_epi32(_mm_sub_epi32(_mm_and_si128(A8, mask4), _mm_and_si128(B8, mask4)));
        diff4 = _mm_add_epi32(diff4, _mm_srai_epi32(_mm_mul_epi32(d2, d2), 12));
    }
    __int64 retval = diff[0];
    retval += diff[1];
    retval += diff[2];
    retval += diff[3];
    return retval;
}


// -----------------------------------------------------------
// Application initialization
// Load a previously saved generation, if available.
// -----------------------------------------------------------
void Game::Init()
{
    for (int i = 0; i < LINES; i++)
        MutateLine(i);
    FILE* f = fopen(LINEFILE, "rb");
    if (f)
    {
        fread(lx1, 4, LINES, f);
        fread(ly1, 4, LINES, f);
        fread(lx2, 4, LINES, f);
        fread(ly2, 4, LINES, f);
        fclose(f);
    }
    Surface* reference = new Surface("assets/image3.png");
    backup = new Surface(SCRWIDTH, SCRHEIGHT);
    ref8 = (BYTE*)MALLOC64(dims);
    for (int i = 0; i < dims; i++)
        ref8[i] = (BYTE)(reference->pixels[i] & 255);

    col8 = (BYTE*)MALLOC64(dims);
    back8 = (BYTE*)MALLOC64(dims);


    fitness = 512 * 512 * 16;
}

// -----------------------------------------------------------
// Application termination
// Save the current generation, so we can continue later.
// -----------------------------------------------------------
void Game::Shutdown()
{
    FILE* f = fopen(LINEFILE, "wb");
    fwrite(lx1, 4, LINES, f);
    fwrite(ly1, 4, LINES, f);
    fwrite(lx2, 4, LINES, f);
    fwrite(ly2, 4, LINES, f);
    fclose(f);
}

// -----------------------------------------------------------
// Main application tick function
// -----------------------------------------------------------
void Game::Tick(float _DT)
{
    tm.reset();
    int lineCount = 0;
    // draw up to lidx
    //memset(screen->pixels, 255, dims_mult);
    memset(col8, 255, dims);
    for (int j = 0; j < lidx; j++, lineCount++)
    {
        unsigned int c = (j << 7) >> LINES_SHIFT;
        DrawWuLine(col8, lx1[j], ly1[j], lx2[j], ly2[j], (c));
    }
    int base = lidx;
    //screen->CopyTo(backup, 0, 0);
    //screen->CopyTo(backup, 0, 0);
    memcpy(back8, col8, dims);

    // iterate and draw from lidx to end
    for (int k = 0; k < ITERATIONS; k++)
    {
        //memcpy(screen->pixels, backup->pixels, dims_mult);
        memcpy(col8, back8, dims);

        MutateLine(lidx);
        for (int j = base; j < LINES; j++, lineCount++)
        {
            unsigned int c = (j << 7) >> LINES_SHIFT;
            DrawWuLine(col8, lx1[j], ly1[j], lx2[j], ly2[j], (c));
        }


        __int64 diff = Evaluate();
        if (diff < fitness)
            fitness = diff;
        else
            UndoMutation(lidx);
        lidx = (lidx + 1) % LINES;
    }


    for (int pos = 0; pos < dims; pos++) {
        uint8_t col = col8[pos];
        screen->pixels[pos] = RGB(col, col, col);
    }



    // stats
    char t[128];
    float elapsed = tm.elapsed();
    float lps = (float)lineCount / elapsed;
    peak = max(lps, peak);
    sprintf(t, "fitness: %i", fitness);
    screen->Bar(0, height_33, 130, height_1, 0);
    screen->Print(t, 2, height_24, 0xffffff);
    sprintf(t, "lps:     %5.2fK", lps);
    screen->Print(t, 2, height_16, 0xffffff);
    sprintf(t, "ips:     %5.2f", iterations_1000 / elapsed);
    screen->Print(t, 2, height_8, 0xffffff);
    sprintf(t, "peak:    %5.2f", peak);
    screen->Print(t, 2, SCRHEIGHT - 32, 0xffffff);
}