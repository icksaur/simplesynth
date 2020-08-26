#include "gui.h"

static void clamp(float & v, float min, float max) { v = v < min ? min : (v > max ? max : v); }
static void clamp(int & v, int min, int max) { v = v < min ? min : (v > max ? max : v); }

bool Rect::contains(int _x, int _y) {
    int hw = w / 2;
    if (abs(x + hw - _x) > hw) return false;
    int hh = h / 2;
    if (abs(y + hh - _y) > hh) return false;
    return true;
}

GUI::GUI() : f() {
    elementColor = 0xDDDDDDFF;
    elementHighlightColor = 0xFFFFFFFF;
    elementBorderColor = sliderHandleColor = textColor = 0x000000FF;
    hotColor = 0xFF5555FF;
    windowColor = 0xBBBBFFFF;
    sliderWidth = 5;
    minDragSize = 10;
    dragThreshold = 2;
    windowBarHeight = 20;
    closed = false;
}
void GUI::pushMouseDown(int which, int x, int y) {
    if (!events.empty() && (events.back().type == mousedown || events.back().type == motion)) {
        Event & c = events.back(); // pitch clicks and motion
        c.type = mousedown;
        c.rect.x = x;
        c.rect.y = y;
    } else {
        events.push_back(Event{ mousedown, (char)which, x, y });
    }
}
void GUI::pushMouseUp(int which, int x, int y) {
    if (which != 0) return; 
    for (int i = events.size() - 1; i >= 0; i--) {
        if (events[i].type == drag) { // end drag
            events.erase(events.begin()+i); // safe because we're iterating backwards
            return;
        } else if (events[i].type == mousedown) { // click
            events[i].type = click; // change mousedown in stream 
            events[i].rect.w = events[i].rect.x;
            events[i].rect.h = events[i].rect.y;
            events[i].rect.x = x;
            events[i].rect.y = y;
        }
    }
}
void GUI::pushMouseMotion(int x, int y) {
    if (events.size()) {
        GUI::Event & e = events.back();
        if (e.type != motion) events.push_back(GUI::Event{ motion, 0, Rect{x, y} });
        else {
            e.rect.x = x;
            e.rect.y = y;
        }
    } else events.push_back(GUI::Event{ motion, 0, Rect{x, y} });
}
void GUI::alignText(int & w, int & h, int & tx, int & ty, int x, int y, const char * label) {
    int tw = f.textWidth(label);
    int th = f.textHeight();
    tx = x + w / 2 - tw / 2; // center horizontal
    ty = y + h / 2 - th / 2; // center vertical
}
void GUI::makeScoped(int & x, int & y) {
    x = scope.x + x;
    y = scope.y + y;
}
// only handles last click in input stream
bool GUI::button(int x, int y, int w, int h, const char * label) {
    makeScoped(x, y);
    int tx, ty;
    alignText(w, h, tx, ty, x, y, label);
    Rect r{ x, y, w, h };
    bool hover = false;
    bool clicked = false;
    int mx, my;
    if (tryGetMouse(mx, my) && r.contains(mx, my)) {
        for (size_t i = 0; i < events.size(); i++) {
            Event & click = events[i];
            if (click.type != GUI::click) continue;
            if (r.contains(click.rect.x, click.rect.y) &&
                r.contains(click.rect.w, click.rect.h)) {
                events.erase(events.begin() + i); // consume the click
                clicked = true;
            }
        }
        hover = r.contains(mx, my);
    }
    f.setColor(elementBorderColor).drawRect(x - 1, y - 1, w + 2, h + 2); // border
    unsigned buttonColor = clicked ? hotColor : hover ? elementHighlightColor : elementColor;
    f.setColor(buttonColor).drawRect(x, y, w, h); // button 
    f.setColor(textColor).drawText(tx, ty, label); // label
    return clicked;
}
void GUI::label(int x, int y, const char * text) {
    makeScoped(x, y);
    f.setColor(textColor).drawText(x, y - f.textHeight()/2, text);
}
bool GUI::checkbox(int x, int y, int w, int h, bool & v) {
    bool toggled = button(x, y, w, h, v ? checkBoxString : "");
    if (toggled) v = !v;
    return toggled;
}
bool GUI::tryGetMouse(int & x, int & y) {
    if (events.empty()) return false;
    x = events.back().rect.x;
    y = events.back().rect.y;
    return true;
}

bool GUI::slider(float & v, int x, int y, int w, int h) {
    makeScoped(x, y);
    if (w < sliderWidth * 2) w = sliderWidth * 2;
    Rect r{ x, y, w, h };
    f.setColor(elementBorderColor).drawRect(x - 1, y - 1, w + 2, h + 2); // border
    bool hot = false;
    bool hover = false;
    int mx, my;
    if (tryGetMouse(mx, my)) {
        hover = r.contains(mx, my);
        GUI::Event & e = events.front();
        if (e.type == mousedown
            && r.contains(e.rect.x, e.rect.y)) {
            v = float(mx - x) / (float)w;
            hot = true;
        }
        else if (e.type == click
            && r.contains(e.rect.w, e.rect.h)) {
            v = float(e.rect.x - x) / (float)w;
            events.erase(events.begin());
            hot = true;
        }
    }
    clamp(v, 0, 1);
    f.setColor(hot ? hotColor : hover ? elementHighlightColor : elementColor);
    f.drawRect(x, y, w, h);
    int valueRange = w - sliderWidth;
    int sliderValue = x + valueRange * v;
    f.setColor(sliderHandleColor);
    f.drawRect(sliderValue, y, sliderWidth, h);
    return hot;
}
bool GUI::dragBox(int & x, int & y, int w, int h) {
    bool hot = false;
    Rect r{ x, y, w, h };
    int mx, my;
    bool hover = false;
    if (tryGetMouse(mx, my)) {
        GUI::Event & e = events.front();
        if (e.type == mousedown) {
            if (r.contains(e.rect.x, e.rect.y)
                && (abs(mx - e.rect.x) > dragThreshold
                    || abs(my - e.rect.y) > dragThreshold)) { // drag started
                e.type = GUI::drag;
                e.rect.w = e.rect.x - x;
                e.rect.h = e.rect.y - y;
                x = mx - e.rect.w;
                y = my - e.rect.h;
                hot = true;
            }
        } else if (e.type == drag) { // drag found
            e.rect.x = mx;
            e.rect.y = my;
            x = mx - e.rect.w;
            y = my - e.rect.h;
            hot = true;
        }
        hover = r.contains(mx, my);
    }
    f.setColor(elementBorderColor).drawRect(x - 1, y - 1, w + 2, h + 2);
    f.setColor(hot ? hotColor : hover ? elementHighlightColor : elementColor).drawRect(x, y, w, h); // box
    return hot;
}
void GUI::line(int x0, int y0, int x1, int y1) {
    int w = x1 - x0;
    int h = y1 - y0;
    if (abs(w) >= abs(h)) { // horizontal box sheared up
        f.setColor(elementHighlightColor).drawRect(x0, y0, w, 1, 0.0f, (float)h/(float)w);
    } else {
        f.setColor(elementHighlightColor).drawRect(x0, y0, 1, h, (float)w/(float)h, 0.0f);
    }
}
void GUI::start(int & x, int & y, int w, int h, const char * text) {
    scope = Rect{ x, y, w, h };
    int dragy = y + h - windowBarHeight;
    dragBox(x, dragy, w - windowBarHeight, windowBarHeight); // draggable title bar
    scope = Rect{ x, dragy - h + windowBarHeight, w, h};
    f.setColor(elementBorderColor).drawRect(scope.x-1, scope.y-1, w+2, h - windowBarHeight); // window border
    closed = button(
        w - windowBarHeight + 1, h - windowBarHeight,
        windowBarHeight - 1, windowBarHeight, "X");
    y = dragy - h + windowBarHeight;
    label(0, h - windowBarHeight/2, text); // title
    f.setColor(windowColor).drawRect(x, y, w, h - windowBarHeight - 1); // window client area
}
bool GUI::end() {
    bool cleaned = false;
    for (auto i = events.begin(); i != events.end(); i++) {
        if (i->type != click) {
            events.erase(events.begin(), i);
            return closed;
        }
    }
    events.clear(); // nothing but unconsumed clicks
    return closed;
}
void GUI::render() {
    f.render();
    f.reset(); // remove this to render many times without more input
}
