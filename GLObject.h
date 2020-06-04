#pragma once

#include <vector>
#include <string>

void glThrowOnError();

// create a valid shader program given a combined vertex+fragment source.
// vertex source must be wrapped with #ifdef VERTEX_SHADER
// frament source must be wrapped with #ifdef FRAGMENT_SHADER
class ShaderProgram {
    unsigned program;
public:
    ShaderProgram(const char * sources, unsigned length);
    ShaderProgram(const char * sources);
    ~ShaderProgram();
    void use();
    int getUniformLocation(const char * name);
    int getAttribLocation(const char * name);
    unsigned getProgram() const;
};

class Texture {
    unsigned dimension;
    unsigned texture;
public:
    Texture(unsigned dimension);
    ~Texture();
    void bind();
    unsigned getDimensions() const;
};

class ArrayTexture {
    unsigned dimension, depth;
    unsigned texture;
public:
    ArrayTexture(unsigned dimension, unsigned depth, const std::vector<char> & bytes);
    ArrayTexture(unsigned dimension, unsigned depth, void * bytes, size_t length);
    ~ArrayTexture();
    void bind();
    unsigned getDimensions() const;
    unsigned getDepth() const;
};

struct AttributeDefinition {
    int index;
    unsigned short size;
    int type;
    int divisor;
};

class VAOBuilder {
    std::vector<std::pair<unsigned short, AttributeDefinition> > attributes;
    unsigned short vboCount;
    int divisor;
public:
    VAOBuilder(unsigned short vboCount);
    VAOBuilder & addFloats(unsigned short vboIndex, int attributeIndex, unsigned short size);
    VAOBuilder & addInts(unsigned short vboIndex, int attributeIndex, unsigned short size);
    VAOBuilder & addUInts(unsigned short vboIndex, int attributeIndex, unsigned short size);
    void setDivisor(int divisor);
    void setAttributes(unsigned short vboIndex, class ShaderProgram & program);
    unsigned short getVboCount() const;
};

class VAO {
    unsigned vao, ibo;
    std::vector<unsigned> vbos;
    ShaderProgram & program;
    VAO() = delete;
public:
    VAO(class ShaderProgram & program, VAOBuilder & builder);
    void vertexData(void * data, unsigned size, unsigned whichBuffer);
    void vertexData(void * data, unsigned size);
    void indexData(void * data, unsigned size);
    void bind();
    void unbind();
    ShaderProgram * getProgram() const;
};

class FBO {
    unsigned fbo, rbo;
    Texture texture;
public:
    FBO(unsigned dimension);
    ~FBO();
    void bind();
    void bindTexture();
    static void unbind();
};

struct UBOCollection {
    unsigned next;
    unsigned short bytes;
    std::vector<unsigned> ubos;
};

class Uncopyable {
public:
    Uncopyable() = default;
    Uncopyable(Uncopyable&) = delete;
    Uncopyable(Uncopyable&&) = delete;
    Uncopyable & operator=(Uncopyable&) = delete;
};

class UBO : public Uncopyable {
    size_t size;
    unsigned object;
public:
    UBO(void * data, size_t size);
    void setData(void * data, size_t size);
    void bind(int bindingPoint);
    ~UBO();
};

class UBOProvider {
    int sizeToCollection[sizeof(unsigned short)];
    std::vector<UBOCollection*> collections;
public: 
    UBOProvider();
    UBOProvider(const UBOProvider &) = delete;
    UBOProvider & operator=(const UBOProvider &) = delete;
    ~UBOProvider();

    static UBOProvider & instance();

    unsigned getUBO(void * bytes, unsigned short size);
    void clear();
};
