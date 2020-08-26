#pragma once

#include "BitFont.h"

struct Rect {
    int x, y, w, h;
    bool contains(int _x, int _y);
};

class GUI {
public:
    GUI::GUI();
    void pushMouseDown(int which, int x, int y);
    void pushMouseUp(int which, int x, int y);
    void pushMouseMotion(int x, int y);
    void start(int & x, int & y, int w, int h, const char * text);
    bool end();
    bool checkbox(int x, int y, int w, int h, bool & v);
    bool button(int x, int y, int w, int h, const char * label);
    void label(int x, int y, const char * text);
    void line(int x0, int y0, int x1, int y1);
    bool tryGetMouse(int & x, int & y);
    bool slider(float & v, int x, int y, int w, int h);
    void render();
private:
    int elementColor, elementBorderColor, sliderHandleColor, elementHighlightColor,
        hotColor, textColor, windowColor;
    int windowBarHeight;
    int sliderWidth;
    int minDragSize;
    int dragThreshold;
    const char * checkBoxString = "x";
    enum EventType {
        mousedown, click, drag, motion
    };
    struct Event {
        EventType type;
        char which;
        Rect rect;
    };
    BitFont f;
    std::vector<Event> events;
    GUI::Rect scope;
    bool closed;

    void alignText(int & w, int & h, int & tx, int & ty, int x, int y, const char * label);
    void makeScoped(int & x, int & y);
    bool dragBox(int & x, int & y, int w, int h);
};
