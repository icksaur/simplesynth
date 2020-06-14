#include "GLObject.h"
#include "glxw.h"

#include <vector>

void buildShaderProgram(unsigned & shaderProgram, const char * source, unsigned length);

ShaderProgram::ShaderProgram(const char * source) {
    buildShaderProgram(program, source, strlen(source)); // throws on build failure
}

ShaderProgram::ShaderProgram(const char * source, unsigned length) {
    buildShaderProgram(program, source, length); // throws on build failure
}

void ShaderProgram::use() {
    glUseProgram(program);
}

int ShaderProgram::getUniformLocation(const char * name) {
    return glGetUniformLocation(program, name);
}

int ShaderProgram::getAttribLocation(const char * name) {
    return glGetAttribLocation(program, name);
}

unsigned ShaderProgram::getProgram() const {
    return program;
}

ShaderProgram::~ShaderProgram() {
    glDeleteProgram(program);
}

ArrayTexture::ArrayTexture(unsigned dimension, unsigned depth, void * bytes, size_t size)
    : dimension(dimension), depth(depth)
{
    if (dimension == 0 || depth == 0) throw std::runtime_error("texture size is zero");
    unsigned levelSize = dimension * dimension * 4;
    if (levelSize * depth != size) throw std::runtime_error("bytes are not 4x pixel count");

    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D_ARRAY, texture);
    glTexStorage3D(GL_TEXTURE_2D_ARRAY, 1, GL_RGBA8, dimension, dimension, depth);

    // reverse image order
    for (unsigned i = 0; i < depth; i++) {
        glTexSubImage3D(
            GL_TEXTURE_2D_ARRAY,
            0, // level
            0, 0, // x and y offsets
            depth - i - 1, // z offset
            dimension, dimension,
            1, // depth
            GL_BGRA, GL_UNSIGNED_BYTE,
            &((char*)bytes)[i * levelSize]); // skip 
    }

    // common defaults
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glThrowOnError();
}

ArrayTexture::ArrayTexture(unsigned dimension, unsigned depth, const std::vector<char> & bytes)
    : ArrayTexture(dimension, depth, (void*)&bytes[0], bytes.size()) {
}

ArrayTexture::~ArrayTexture() {
    glDeleteTextures(1, &texture);
}

void ArrayTexture::bind() {
    glBindTexture(GL_TEXTURE_2D_ARRAY, texture);
}

unsigned ArrayTexture::getDimensions() const {
    return dimension;
}

unsigned ArrayTexture::getDepth() const {
    return depth;
}

Texture::Texture(unsigned dimension) {
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, dimension, dimension, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
    glThrowOnError();
}

Texture::~Texture() {
    glDeleteTextures(1, &texture);
}

void Texture::bind() {
    glBindTexture(GL_TEXTURE_2D, texture);
}

unsigned Texture::getDimensions() const {
    return dimension;
}

/*
const char * GlobalShaderHeader = R"(
#version 420 core
layout(binding=15, std140) uniform global_data {
    float time;
    float deltaTime;
    ivec2 resolution;
};   
layout(binding=14, std140) uniform batch_data {
    mat4 viewProjection;
};
)";
*/
const char * GlobalShaderHeader = R"(#version 430 core
)";

const char * VertexShaderHeader = R"(#define VERTEX_SHADER 1
)";

const char * FragmentShaderHeader = R"(#define FRAGMENT_SHADER 1
)";

unsigned compileShaderObject(unsigned & object, GLenum shaderType, const char * source, unsigned length) {
    const char * header = nullptr;
    if (shaderType == GL_VERTEX_SHADER) {
        header = VertexShaderHeader;
    }
    else if (shaderType == GL_FRAGMENT_SHADER) {
        header = FragmentShaderHeader;
    }

    glThrowOnError();

    const char * sources[3]{ GlobalShaderHeader, header, source };
    int lengths[3]{ (int)strlen(GlobalShaderHeader), (int)strlen(header), (int)length };
    glShaderSource(object, 3, sources, lengths);
    glCompileShader(object);

    int compileStatus;
    glGetShaderiv(object, GL_COMPILE_STATUS, &compileStatus);

    if (!compileStatus) { // failure
        int bufferSize, logLength;
        glGetShaderiv(object, GL_INFO_LOG_LENGTH, &bufferSize);
        std::string compileError;
        compileError.resize(bufferSize);
        glGetShaderInfoLog(object, bufferSize, &logLength, &compileError[0]);
        glDeleteShader(object);
        throw std::runtime_error(compileError);
    }

    return object;
}

void buildShaderProgram(unsigned & shaderProgram, const char * source, unsigned length) {
    unsigned vertexShader;
    unsigned fragmentShader;
    if (!glIsProgram(shaderProgram)) {
        shaderProgram = glCreateProgram();
        vertexShader = glCreateShader(GL_VERTEX_SHADER);
        fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glAttachShader(shaderProgram, vertexShader);
        glAttachShader(shaderProgram, fragmentShader);
    }
    else {
        int count;
        unsigned shaders[2];
        glGetAttachedShaders(shaderProgram, 2, &count, shaders);
        int type;
        for (int i = 0; i < 2; i++) {
            glGetShaderiv(shaders[i], GL_SHADER_TYPE, &type);
            switch (type) {
            case GL_VERTEX_SHADER: vertexShader = shaders[i]; break;
            case GL_FRAGMENT_SHADER: fragmentShader = shaders[i]; break;
            default: break;
            }
        }
    }

    compileShaderObject(vertexShader, GL_VERTEX_SHADER, source, length);
    compileShaderObject(fragmentShader, GL_FRAGMENT_SHADER, source, length);
    glThrowOnError();

    glLinkProgram(shaderProgram);
    glThrowOnError();
    glDeleteShader(vertexShader);
    glThrowOnError();
    glDeleteShader(fragmentShader);
    glThrowOnError();

    int linkStatus;
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &linkStatus);

    if (!linkStatus) { // failure
        int logLength, bufferSize;
        std::string linkError;
        glGetProgramiv(shaderProgram, GL_INFO_LOG_LENGTH, &bufferSize);
        linkError.resize(bufferSize);
        glGetProgramInfoLog(shaderProgram, bufferSize, &logLength, &linkError[0]);
        throw std::runtime_error(linkError);
    }
}

const int MAXVBOS = 2;

VAOBuilder::VAOBuilder(unsigned short vboCount) : divisor(0), vboCount(vboCount) {
    if (vboCount == 0) throw std::runtime_error("VAO requires at least 1 VBO");
    if (vboCount > MAXVBOS) throw std::runtime_error("too many VBOs for a VAO");
}

VAOBuilder & VAOBuilder::addFloats(unsigned short vboIndex, int attributeIndex, unsigned short size) {
    if (vboIndex >= vboCount) throw std::runtime_error("VBO index higher than intended count");
    attributes.push_back(std::make_pair(vboIndex, AttributeDefinition{ attributeIndex, size, GL_FLOAT, divisor }));
    return *this;
}
VAOBuilder & VAOBuilder::addInts(unsigned short vboIndex, int attributeIndex, unsigned short size) {
    if (vboIndex >= vboCount) throw std::runtime_error("VBO index higher than intended count");
    attributes.push_back(std::make_pair(vboIndex, AttributeDefinition{ attributeIndex, size, GL_INT, divisor }));
    return *this;
}
VAOBuilder & VAOBuilder::addUInts(unsigned short vboIndex, int attributeIndex, unsigned short size) {
    if (vboIndex >= vboCount) throw std::runtime_error("VBO index higher than intended count");
    attributes.push_back(std::make_pair(vboIndex, AttributeDefinition{ attributeIndex, size, GL_UNSIGNED_INT, divisor }));
    return *this;
}
void VAOBuilder::setDivisor(int divisor) {
    this->divisor = divisor;
}

void VAOBuilder::setAttributes(unsigned short vboIndex, ShaderProgram & program) {
    program.use();
    unsigned short vboVertexSize = 0; // aka stride
    for (std::pair<unsigned short, AttributeDefinition> & attributePair : attributes) {
        if (attributePair.first == vboIndex) {
            vboVertexSize += attributePair.second.size * (attributePair.second.type == GL_INT ? sizeof(int) : sizeof(float));
        }
    }

    short offset = 0;
    for (std::pair<unsigned short, AttributeDefinition> & attributePair : attributes) {
        if (attributePair.first != vboIndex) {
            continue; // only 
        }
        int location = attributePair.second.index;
        if (location == -1) throw std::runtime_error(std::string("attribute location is invalid"));
        glEnableVertexAttribArray(location);
        if (attributePair.second.type == GL_FLOAT) {
            glVertexAttribPointer(
                location,
                attributePair.second.size,
                attributePair.second.type,
                GL_FALSE,
                vboVertexSize,
                (void*)offset);
        } else {
            // This is NOT the same as AttribPointer with GL_FALSE. The Above turns integers into floats.
            glVertexAttribIPointer(
                location,
                attributePair.second.size,
                attributePair.second.type,
                vboVertexSize,
                (void*)offset);
        }
        glVertexAttribDivisor(location, attributePair.second.divisor); // 0 means per vertex, n means per attribute
        offset += attributePair.second.size * (attributePair.second.type == GL_INT ? sizeof(int) : sizeof(float));
    }
}

unsigned short VAOBuilder::getVboCount() const {
    return vboCount;
}

VAO::VAO(ShaderProgram & program, VAOBuilder & builder) : program(program) {
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &ibo);
    vbos.resize(builder.getVboCount());
    glGenBuffers(builder.getVboCount(), &vbos[0]);

    struct CleanUp {
        unsigned vao, ibo, *vbos, vboCount;
        bool success;
        ~CleanUp() {
            if (success) return;
            glDeleteVertexArrays(1, &vao);
            glDeleteBuffers(vboCount, vbos);
            glDeleteBuffers(1, &ibo);
        }
    } cleanup = { vao, ibo, &vbos[0], builder.getVboCount(), false };

    program.use();
    glBindVertexArray(vao);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    for (unsigned short i = 0; i < builder.getVboCount(); i++) {
        glBindBuffer(GL_ARRAY_BUFFER, vbos[i]);
        builder.setAttributes(i, program);
    }
    glBindVertexArray(0);

    cleanup.success = true;
}

void VAO::vertexData(void * data, unsigned size) {
    vertexData(data, size, 0);
}

void VAO::vertexData(void * data, unsigned size, unsigned whichBuffer) {
    if (whichBuffer >= vbos.size()) throw std::runtime_error("no such VBO in this VAO");
    glBindBuffer(GL_ARRAY_BUFFER, vbos[whichBuffer]);
    glBufferData(GL_ARRAY_BUFFER, size, data, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void VAO::indexData(void * data, unsigned size) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, size, data, GL_STATIC_DRAW);
}

void VAO::bind() {
    glBindVertexArray(vao);
}

void VAO::unbind() {
    glBindVertexArray(0);
}

ShaderProgram * VAO::getProgram() const {
    return &program;
}

FBO::FBO(unsigned dimension) :texture(dimension) {
    glGenRenderbuffers(1, &rbo);
    glBindRenderbuffer(GL_RENDERBUFFER, rbo);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT16, dimension, dimension);
    glBindRenderbuffer(GL_RENDERBUFFER, 0); // maybe not needed

    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, fbo, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rbo);
    glThrowOnError();
}

FBO::~FBO() {
    glDeleteRenderbuffers(1, &rbo);
    glDeleteFramebuffers(1, &fbo);
}

void FBO::bind() {
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
}

void FBO::bindTexture() {
    texture.bind();
}

void FBO::unbind() {
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

UBOProvider::UBOProvider() {
    memset(sizeToCollection, -1, 256);
}

UBOProvider::~UBOProvider() {
    for (UBOCollection * collection : collections) {
        glDeleteBuffers(collection->ubos.size(), &(collection->ubos[0]));
    }
}

UBOProvider & UBOProvider::instance() {
    static UBOProvider provider;
    return provider;
}

//SETTING(maxUboSize, 5000);
const int maxUboSize = 5000;

unsigned UBOProvider::getUBO(void * bytes, unsigned short size) {
    if (size > maxUboSize) throw std::runtime_error("max ubo size exceeded");
    int index = sizeToCollection[size];
    if (index < 0) { // no collection for this size.  Make it.
        UBOCollection * newCollection = new UBOCollection;
        collections.push_back(newCollection);
        newCollection->next = -1;
        newCollection->bytes = size;
        sizeToCollection[size] = collections.size()-1;
        index = sizeToCollection[size];
    }

    if (index < 0) throw std::runtime_error("ubo size collection not found");

    unsigned ubo;
    UBOCollection & coll = *collections[index];
    if (coll.next < coll.ubos.size()) {
        ubo = coll.ubos[coll.next++];
    } else {
        // make a new ubo of right size
        glGenBuffers(1, &ubo);
        coll.ubos.push_back(ubo);
        coll.next++;
    }
    glBindBuffer(GL_UNIFORM_BUFFER, ubo);
    glBufferData(GL_UNIFORM_BUFFER, size, bytes, GL_DYNAMIC_DRAW);
    return ubo;
}

void UBOProvider::clear() {
    for (UBOCollection * collection : collections) {
        collection->next = 0;
    }
}

void setUboData(unsigned object, void * data, size_t size) {
    glBindBuffer(GL_UNIFORM_BUFFER, object);
    glBufferData(GL_UNIFORM_BUFFER, size, data, GL_DYNAMIC_DRAW);
}

UBO::UBO(void * data, size_t size):size(size) {
    if ((int)size > maxUboSize) throw std::runtime_error("max ubo size exceeded");
    glGenBuffers(1, &object);
    setUboData(object, data, size);
}

UBO::~UBO(){
    glDeleteBuffers(1, &object);
}

void UBO::setData(void * data, size_t size) {
    if (size != this->size) throw std::runtime_error("UBO size changed");
    setUboData(object, data, size);
}

void UBO::bind(int bindingPoint) {
    glBindBufferBase(GL_UNIFORM_BUFFER, bindingPoint, object);
}

#define ERRORCASE(CODE) case CODE: return #CODE;

static const char * getError(GLuint error) {
    switch (error) {
        ERRORCASE(GL_NO_ERROR);
        ERRORCASE(GL_INVALID_ENUM);
        ERRORCASE(GL_INVALID_VALUE);
        ERRORCASE(GL_INVALID_OPERATION);
        ERRORCASE(GL_OUT_OF_MEMORY);
        ERRORCASE(GL_INVALID_FRAMEBUFFER_OPERATION);
        ERRORCASE(GL_STACK_OVERFLOW);
        ERRORCASE(GL_STACK_UNDERFLOW);
    }
    return "unknown error";
}

void glThrowOnError() {
#if 1
    GLenum error = glGetError();
    if (GL_NO_ERROR != error) {
        printf("0x%x: %s", error, getError(error));
#ifdef _DEBUG
        DebugBreak();
#else
        throw std::runtime_error("OpenGL Error");
#endif
    }
#endif
}
