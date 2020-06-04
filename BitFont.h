#pragma once

#include <vector>

struct BitFont {
    struct Glyph {
        float x, y, w, h;
    };
    struct Instance {
        float x, y;
        float w, h;
        float r, g, b, a;
        int glyphId;
    };
    float r, g, b, a;
    class ArrayTexture * texture;
    class ShaderProgram * program;
    class VAO * vao;
    class UBO * ubo;
    unsigned char asciiToGlyph[256];
    bool dirty;
    std::vector<BitFont::Glyph> glyphs;
    std::vector<BitFont::Instance> instances;
    const int scale = 2;
    int xrez, yrez;

    BitFont(); 
    ~BitFont();
    BitFont & drawRect(int x, int y, int w, int h);
    BitFont & drawText(int x, int y, const char * text);
    void patchRect(int previous, int x, int y);
    int textWidth(const char * text);
    int textHeight();
    void reset();
    BitFont & setColor(float _r, float _g, float _b, float _a);
    BitFont & setColor(unsigned color);
    void setResolution(int xrez, int yrez);
    void render(int first = 0, int count = -1);
};
