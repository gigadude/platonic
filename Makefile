TARGET := platonic

LIBS := -L/usr/X11R6/lib -lglut -lGLU -lGL -lXmu -lXext -lX11 -lm
OPTS := -O6 -ffast-math -mfpmath=sse -msse2
#OPTS := -g -ffast-math -mfpmath=sse -msse2

all: $(TARGET)

%: %.cpp
	g++ -o $@ $(OPTS) $< $(LIBS)

%: %.c
	g++ -o $@ $(OPTS) $< $(LIBS)

%.pm: %.yp
	eyapp -l -s -v $<

clean: $(wildcard $(TARGET) glslparser.pm glslparser.output shaders.h *.o)
	rm -f $^

$(TARGET): ssevector.h glslglue.h shaders.h

shaders.h: glueglsl glslparser.pm shaders.glsl
	glueglsl shaders.glsl > $@
