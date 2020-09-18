#include "GLObject.h"
#include "BitFont.h"

#include <GLXW/glxw.h>

#include <stdexcept>

// location in a texture of a glyph where 0,0 is the upper-left corner of the texture
struct Glyph {
    int x; // left side of glyph
    int y; // top of glyph
    int w; // width in pixels
    int h; // height in pixels
    int xoffset; // amount to move glyph right (unused?)
    int yoffset; // descender height
    int xadvance; // how far forward for the next glyph
};

const int ascenderHeight = 8;
const int descenderHeight = 2;
const int fontHeight = ascenderHeight + descenderHeight;

// x, y, w, h, xoff, yoff, advance
Glyph bitmapFontGlyphs[95] = {
    { 32, 37, 6, 6, 0, 0, 5 }, // a
    { 57, 4, 6, 7, 0, 0, 5 },
    { 30, 57, 5, 6, 0, 0, 4 },
    { 48, 8, 6, 7, 0, 0, 5 },
    { 54, 11, 6, 6, 0, 0, 5 },
    { 44, 15, 5, 7, 0, 0, 4 },
    { 49, 17, 6, 8, 0, 2, 5 },
    { 43, 22, 6, 7, 0, 0, 5 },
    { 27, 22, 3, 7, 0, 0, 2 },
    { 55, 17, 4, 9, 0, 2, 3 },
    { 49, 25, 6, 7, 0, 0, 5 },
    { 60, 11, 3, 7, 0, 0, 2 },
    { 55, 26, 7, 6, 0, 0, 6 },
    { 38, 29, 6, 6, 0, 0, 5 },
    { 44, 32, 6, 6, 0, 0, 5 },
    { 38, 35, 6, 8, 0, 2, 5 },
    { 50, 32, 6, 8, 0, 2, 5 },
    { 44, 38, 5, 6, 0, 0, 4 },
    { 56, 32, 6, 6, 0, 0, 5 },
    { 56, 38, 5, 7, 0, 0, 4 },
    { 49, 40, 6, 6, 0, 0, 5 },
    { 55, 45, 6, 6, 0, 0, 5 },
    { 35, 43, 7, 6, 0, 0, 6 },
    { 42, 44, 5, 6, 0, 0, 4 },
    { 35, 49, 6, 8, 0, 2, 5 },
    { 35, 57, 6, 6, 0, 0, 5 }, // z
    { 19, 34, 6, 7, 0, 0, 5 }, // A
    { 20, 41, 6, 7, 0, 0, 5 },
    { 24, 48, 5, 7, 0, 0, 4 },
    { 24, 55, 6, 7, 0, 0, 5 },
    { 20, 1, 5, 7, 0, 0, 4 },
    { 20, 8, 5, 7, 0, 0, 4 },
    { 20, 15, 6, 7, 0, 0, 5 },
    { 21, 22, 6, 7, 0, 0, 5 },
    { 25, 29, 5, 7, 0, 0, 4 },
    { 26, 36, 6, 7, 0, 0, 5 },
    { 29, 43, 6, 7, 0, 0, 5 },
    { 30, 50, 5, 7, 0, 0, 4 },
    { 25, 1, 7, 7, 0, 0, 6 },
    { 25, 8, 6, 7, 0, 0, 5 },
    { 26, 15, 6, 7, 0, 0, 5 },
    { 31, 8, 6, 7, 0, 0, 5 },
    { 30, 22, 6, 8, 0, 1, 5 },
    { 32, 1, 6, 7, 0, 0, 5 },
    { 32, 15, 6, 7, 0, 0, 5 },
    { 37, 8, 5, 7, 0, 0, 4 },
    { 38, 1, 6, 7, 0, 0, 5 },
    { 32, 30, 6, 7, 0, 0, 5 },
    { 36, 22, 7, 7, 0, 0, 6 },
    { 38, 15, 6, 7, 0, 0, 5 },
    { 42, 8, 6, 7, 0, 0, 5 },
    { 44, 1, 5, 7, 0, 0, 4 }, // Z
    { 1, 1, 3, 7, 0, 0, 2 }, // !
    { 17, 43, 3, 6, 0, 1, 2 },
    { 1, 30, 7, 7, 0, 0, 6 },
    { 17, 8, 3, 5, 0, 1, 2 },
    { 14, 20, 6, 7, 0, 0, 5 }, // ?
    { 8, 16, 5, 5, 0, 0, 4 }, // *
    { 8, 9, 4, 7, 0, 0, 3 },
    { 11, 1, 4, 7, 0, 0, 3 },
    { 57, 1, 6, 3, 0, 4, 5 },
    { 12, 8, 5, 5, 0, 1, 4 },
    { 15, 5, 5, 3, 0, 2, 4 },
    { 1, 58, 5, 5, 0, 1, 4 },
    { 4, 5, 3, 3, 0, 0, 2 }, // .
    { 15, 1, 4, 4, 0, 0, 3 },
    { 13, 13, 7, 7, 0, 0, 6 },
    { 152, 46, 3, 7, 0, 0, 2 },
    { 4, 1, 5, 4, 0, 0, 4 },
    { 8, 5, 3, 4, 0, 0, 2 },
    { 14, 27, 7, 7, 0, 0, 6 },
    { 1, 15, 7, 7, 0, 0, 6 },
    { 1, 22, 6, 8, 0, 0, 5 },
    { 6, 58, 5, 4, 0, 0, 4 },
    { 7, 22, 7, 7, 0, 0, 6 },
    { 147, 46, 5, 7, 0, 0, 4 },
    { 141, 50, 5, 7, 0, 0, 4 },
    { 49, 1, 4, 7, 0, 0, 3 },
    { 53, 1, 4, 7, 0, 0, 3 },
    { 1, 8, 7, 7, 0, 0, 6 },
    { 1, 37, 6, 7, 0, 0, 5 },
    { 1, 44, 4, 7, 0, 0, 3 },
    { 1, 51, 6, 7, 0, 0, 5 },
    { 5, 44, 6, 7, 0, 0, 5 },
    { 7, 37, 6, 7, 0, 0, 5 },
    { 8, 29, 6, 7, 0, 0, 5 },
    { 7, 51, 6, 7, 0, 0, 5 },
    { 11, 44, 6, 7, 0, 0, 5 },
    { 13, 36, 6, 7, 0, 0, 5 },
    { 13, 51, 6, 7, 0, 0, 5 },
    { 19, 49, 5, 7, 0, 0, 4 },
    { 19, 56, 5, 7, 0, 0, 4 },
    { 21, 29, 4, 4, 0, 0, 3 },
    { 111, 58, 6, 4, 0, 0, 5 },
    { 0, 0, 3, 0, 0, 0, 4 }, // space
};

int bitmapFontGlyphIndices[95] = {
    97,
    98,
    99,
    100,
    101,
    102,
    103,
    104,
    105,
    106,
    107,
    108,
    109,
    110,
    111,
    112,
    113,
    114,
    115,
    116,
    117,
    118,
    119,
    120,
    121,
    122,
    65,
    66,
    67,
    68,
    69,
    70,
    71,
    72,
    73,
    74,
    75,
    76,
    77,
    78,
    79,
    80,
    81,
    82,
    83,
    84,
    85,
    86,
    87,
    88,
    89,
    90,
    33,
    59,
    37,
    58,
    63,
    42,
    40,
    41,
    95,
    43,
    45,
    61,
    46,
    44,
    47,
    124,
    34,
    39,
    64,
    35,
    36,
    94,
    38,
    123,
    125,
    91,
    93,
    92,
    48,
    49,
    50,
    51,
    52,
    53,
    54,
    55,
    56,
    57,
    60,
    62,
    96,
    126,
    32,
};
//encoded 1bpp image of the Ansi alphabet.  Top row of image is first.
const unsigned imageBits[] = {
    0x0,
    0x0,
    0x0,
    0x0,
    0x25084722,
    0x7127333c,
    0x25048436,
    0x49212100,
    0x2004072a,
    0x49222100,
    0x40422,
    0x71242120,
    0x2448e722,
    0x48c73338,
    0x400000,
    0x24,
    0x0,
    0x24,
    0x20022724,
    0xe3920838,
    0x10270434,
    0x91123800,
    0x842272c,
    0x910e4800,
    0x4400424,
    0xe10248c4,
    0x2400424,
    0x810c3964,
    0x202000,
    0x184,
    0x4000,
    0xc4,
    0x1400838c,
    0x39210004,
    0x3e510412,
    0x41220000,
    0x14220592,
    0x30c71c40,
    0x3e500492,
    0x9222400,
    0x1400038c,
    0x71222440,
    0x1c000,
    0x1c40,
    0x2000,
    0x440,
    0x860c248,
    0xc4481840,
    0x1c800241,
    0x254e0080,
    0x306883c9,
    0x25490000,
    0xc900249,
    0x25492000,
    0x38680248,
    0xc28924f0,
    0x800e000,
    0x200028a8,
    0x11000,
    0x38a8,
    0x797238,
    0x1c024a8,
    0x24415110,
    0x49200000,
    0x470e010,
    0x49200000,
    0x8080010,
    0x51230e38,
    0x10700038,
    0x50049260,
    0x12000600,
    0x20049218,
    0x900,
    0x1c30e70,
    0x18906,
    0x1200200,
    0x18224f02,
    0x39200200,
    0x24618902,
    0x49c50020,
    0x24a24012,
    0x49060070,
    0x24f1800c,
    0x39042420,
    0x18200700,
    0x42420,
    0x480,
    0x2410,
    0x2702,
    0x4a801c00,
    0x338f0482,
    0x8a940000,
    0x10412703,
    0x5080090,
    0x11822002,
    0x85086490,
    0x10440002,
    0x401444a0,
    0x13840030,
    0x8440,
    0x240,
    0x9004400,
    0x441,
    0x9306400,
    0x38618841,
    0x9100000,
    0x4824431,
    0x7080000,
    0x18e1c201,
    0x1100000,
    0x20904001,
    0xc6300000,
    0x3c618070,
    0x0,
    0x848,
    0x0,
    0x448,
    0xcf000000,
    0x38850249,
    0x2000000,
    0x14a0471,
    0x4000000,
    0x38000800,
    0xcf000000,
    0x0,
    0x0,
    0x80000000, // pixel in lower-left for rects
    0x0,
};

BitFont::BitFont() : dirty(true), r(1), g(1), b(1), a(1) {
    int major, minor;
    glGetIntegerv(GL_MAJOR_VERSION, &major);
    glGetIntegerv(GL_MINOR_VERSION, &minor);
    if ((major * 10 + minor) < 43) throw std::runtime_error("OpenGL 4.3 required");

    const char * bitFontShaderSource = R"(
#ifdef VERTEX_SHADER
struct sprite_definition {
vec4 uvco;
};
layout(binding=0, std140) uniform sprite_data {
sprite_definition sprites[];
};
layout(location=0) uniform mat4 viewProjection;
layout(location=0) in vec2 vertex; // object-space vertex of quad, [0,0] to [1,1]
layout(location=1) in vec2 location; // screen-space bottom-left of sprite, from top-left of window at t=0
layout(location=2) in vec2 scale; // screen-space bottom-left of sprite, from top-left of window at t=0
layout(location=3) in vec4 color;
layout(location=4) in vec2 shear; // for drawing lines
layout(location=5) in int index; // sprite in batch
out vec2 textureCoordinate; // 0,0 to 1,1 uv coordiantes from bottomleft of image
out vec4 ocolor;

void main() {
sprite_definition sprite = sprites[index];
ocolor = color;
mat2 shearMatrix = mat2(1.0, shear.y, shear.x, 1.0);
vec2 position = location + (shearMatrix * (scale * vertex));
gl_Position = viewProjection * vec4(position.xy, 0.0, 1.0);
textureCoordinate.x = sprite.uvco.s + sprite.uvco.p * vertex.x;
textureCoordinate.y = sprite.uvco.t + sprite.uvco.q * vertex.y;
}
#endif
#ifdef FRAGMENT_SHADER
uniform sampler2DArray textures;
in vec2 textureCoordinate;
in vec4 ocolor;
void main() {
vec4 texel = texture(textures, vec3(textureCoordinate.st, 0));
gl_FragColor = ocolor * texel;
}
#endif
)";
    program = new ShaderProgram(bitFontShaderSource);
    VAOBuilder builder(2);
    builder.addFloats(0, 0, 2);
    builder.setDivisor(1);
    builder.addFloats(1, 1, 2); // location
    builder.addFloats(1, 2, 2); // scale
    builder.addFloats(1, 3, 4); // color
    builder.addFloats(1, 4, 2); // shear
    builder.addInts(1, 5, 1); // sprite index
    vao = new VAO(*program, builder);
    float square[] = { 0, 0, 1, 0, 1, 1, 0, 1 };
    vao->vertexData(square, sizeof(float) * 8, 0);

    std::vector<unsigned> textureBytes(64 * 64 * 4);

    int npixel = 0;
    for (int row = 63; row >= 0; row--) {
        for (int col = 0; col < 64; col++) {
            int pixel = row * 64 + col;
            unsigned byte4 = imageBits[pixel / 32];

            int bit = 31 - (pixel % 32);

            if ((byte4 >> bit) & 0x1)
                textureBytes[npixel++] = 0xFFFFFFFF; // white
            else
                textureBytes[npixel++] = 0x0;
        }
    }
    texture = new ArrayTexture(64, 1, &textureBytes[0], textureBytes.size());

    BitFont & drawRect(int x, int y, int w, int h);
    const unsigned char questionGlyphIndex = 63;
    memset(asciiToGlyph, questionGlyphIndex, sizeof(asciiToGlyph));
    std::vector<BitFont::Glyph> glyphBytes(96);
    glyphBytes.clear();
    for (int i = 0; i < 95; i++) {
        asciiToGlyph[bitmapFontGlyphIndices[i]] = i;
        ::Glyph & g = bitmapFontGlyphs[i];
        glyphBytes.push_back(BitFont::Glyph{
            (float)g.x / (float)64,
            ((float)64 - (float)g.y - (float)g.h) / (float)64,
            ((float)g.w) / (float)64,
            (float)g.h / (float)64 });
    }
    glyphBytes.push_back(BitFont::Glyph{ 0, 0, (float)1/(float)64, (float)1/(float)64}); // white pixel for rectangles
    ubo = new UBO(&glyphBytes[0], sizeof(BitFont::Glyph) * glyphBytes.size());

    int viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    setResolution(viewport[2], viewport[3]);
}
BitFont::~BitFont() {
    delete program;
    delete texture;
    delete vao;
    delete ubo;
}

BitFont & BitFont::drawRect(int x, int y, int w, int h) {
    return drawRect(x, y, w, h, 0.0f, 0.0f);
}
BitFont & BitFont::drawRect(int x, int y, int w, int h, float sx, float sy) {
    dirty = true;
    instances.push_back(BitFont::Instance{
        (float)(x - xrez/2), (float)(y - yrez/2), // 0,0 to bottomleft
        (float)w, (float)h,
        r, g, b, a,
        sx, sy,
        95 });
    return *this;
}
void BitFont::patchRect(int which, int x, int y) {
    which = instances.size() + which;
    if (which >= 0 && which < (int)instances.size()) {
        Instance & i = instances[which];
        i.x = (float)(x - xrez / 2);
        i.y = (float)(y - yrez / 2);
    }
}
// draw text from 0,0 being bottomleft pixel of given resolution
// y is baseline-descender height, not baseline. (aka bottom of 'p' not 'b')
BitFont & BitFont::drawText(int x, int y, const char * text) {
    dirty = true; 
    int hh = yrez / 2;
    int hw = xrez / 2;
    while (*text != 0) {
        ::Glyph & gl = bitmapFontGlyphs[asciiToGlyph[*text]];
        instances.push_back(BitFont::Instance{
            (float)(x - hw), (float)(y - (gl.yoffset-descenderHeight)*scale - hh), // 0,0 to bottomleft
            (float)gl.w*scale, (float)gl.h*scale,
            r, g, b, a,
            0.0f, 0.0f,
            (int)asciiToGlyph[*text] });
        x += gl.xadvance * (int)scale;
        text++;
    }
    return *this;
}
int BitFont::textWidth(const char * text) {
    int x = 0;
    while (*text != 0) {
        ::Glyph & g = bitmapFontGlyphs[asciiToGlyph[*text]];
        x += g.xadvance*scale;
        text++;
    }
    return x;
}
int BitFont::textHeight() {
    return scale * (ascenderHeight + descenderHeight);
}
void BitFont::reset() {
    dirty = true;
    instances.clear();
}
BitFont & BitFont::setColor(float _r, float _b, float _g, float _a) {
    r = _r; b = _b; g = _g; a = _a;
    return *this;
}
BitFont & BitFont::setColor(unsigned color) {
    unsigned char * c = (unsigned char *)&color;
    r = (float)c[3] / 255.0f;
    g = (float)c[2] / 255.0f;
    b = (float)c[1] / 255.0f;
    a = (float)c[0] / 255.0f;
    return *this;
}
void BitFont::setResolution(int xrez, int yrez) {
    this->xrez = xrez;
    this->yrez = yrez;
    const float hw = (float)xrez / 2, hh = (float)yrez / 2;
    const float znear = -1.0f, zfar = 1.0f;
    float projection[16] = { 0.0f };
    projection[0] = 1.0f / hw;
    projection[5] = 1.0f / hh;
    projection[10] = -2.0f / (zfar - znear);
    projection[14] = -((zfar + znear) / (zfar - znear));
    projection[15] = 1.0f;
    program->use();
    glUniformMatrix4fv(0, 1, GL_TRUE, projection);
}
void BitFont::render(int first, int count) {
    if (instances.size() == 0) return;
    program->use();
    if (dirty) {
        vao->vertexData(&instances[0], sizeof(BitFont::Instance) * instances.size(), 1);
        dirty = false;
    }
    ubo->bind(0);
    glActiveTexture(GL_TEXTURE0 + 0);
    texture->bind();
    vao->bind();
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDrawArraysInstanced(GL_TRIANGLE_FAN, first, 4, count == -1 ? instances.size()-first : count);
}
