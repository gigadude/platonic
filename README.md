Platonic
========

Example OpenGL/GLSL framework, and demo platonic solid renderer.

This little demo shows off an extended GLSL parser (glslparser.yp)
implemented with PERL Parse::Eyapp. The shader extensions allow
you to specify all shaders in one file, using named vertexshader,
geometryshader and fragmentshader blocks. You can also group shaders
into programs with the named program block. The resulting objects
are compiled into a C++ header which encapsulates the interface to
each program.

Note:

All definitions which are outside the block declarations are shared
with all subsequent blocks.

See the shaders.glsl file and the generated shaders.h for an example.

DEPENDENCIES
------------

* Linux
* GL, GLU, freeglut
* Parse::Lex, Parse::Eyapp
The PERL modules can be installed with:

```bash
sudo cpan Parse::Lex Parse::Eyapp
```

TO RUN
------

```bash
make && ./platonic
```

TODO
----

* document the language.
* upgrade to latest GLSL spec.
* the parser gets line #'s in the line directives off-by-one sometimes.
* add support for structures, tesselation shaders, bindless graphics, etc.
* describe platonic a bit.
