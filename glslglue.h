//
// glslglue.h - base classes used by glueglsl
//

#if !defined(__glslglue_h)
#define __glslglue_h

#include <memory>

namespace glsl
{
	static bool strict = false;

	//
	// program base class, generated subclasses contain specific shader and interface instances
	//

	class program
	{
		GLuint _handle;
	public:
		virtual const char *name() const = 0;
		GLuint handle() { return( _handle ); }
		void handle( GLuint h ) { _handle = h; }
		program() { handle( glCreateProgram() ); }
		~program()
		{
			if (handle()) glDeleteProgram( handle() );
		}
		// attach and bind get overridden in generated code
		virtual bool compile( char *errlog, GLsizei *len ) = 0;
		virtual bool init( char *errlog, GLsizei *len ) = 0;
		bool build( char *errlog, GLsizei *len )
		{
			if (!compile( errlog, len ))
				return( false );

			glLinkProgram( handle() );
			GLint st;
			glGetProgramiv( handle(), GL_LINK_STATUS, &st );
			if (!st)
			{
				glGetProgramInfoLog( handle(), *len, len, errlog );
				return( false );
			}

			use();

			if (!init( errlog, len ))
			{
				debug();
				return( false );
			}

			return( true );
		}
		bool info( char *log, GLsizei *len )
		{
			glGetProgramInfoLog( handle(), *len, len, log );
		}
		void debug()
		{
			printf( "program %s:\n", name() );
			GLint attrs;
			glGetProgramiv( handle(), GL_ACTIVE_ATTRIBUTES, &attrs );
			printf( "%d attributes:\n", attrs );
			for (int i = 0; i < attrs; ++i)
			{
				GLint size;
				GLenum type;
				GLchar name[256];
				GLsizei len = sizeof(name);
				glGetActiveAttrib( handle(), i, len, &len, &size, &type, name );
				printf( "%d: %s [%d] 0x%x\n", i, name, size, type );
			}
			GLint uniforms;
			glGetProgramiv( handle(), GL_ACTIVE_UNIFORMS, &uniforms );
			printf( "%d uniforms:\n", uniforms );
			for (int i = 0; i < uniforms; ++i)
			{
				GLint size;
				GLenum type;
				GLchar name[256];
				GLsizei len = sizeof(name);
				glGetActiveUniform( handle(), i, len, &len, &size, &type, name );
				printf( "%d: %s [%d] 0x%x\n", i, name, size, type );
			}
		}
		void use() { glUseProgram( handle() ); }
	};

	//
	// shader base class
	//

	class shader
	{
		const char *text;
		GLuint _handle;

	public:
		virtual const char *name() const = 0;
		GLuint handle() const { return( _handle ); }
		void handle( GLuint h ) { _handle = h; }
		shader( const char *text ) : text(text) { handle( 0 ); }
		~shader()
		{
			if (handle()) glDeleteShader( handle() );
		}
		bool compile( char *errlog, GLsizei *len )
		{
			const GLchar *src[1] = { text };
			glShaderSource( handle(), 1, src, NULL );
			glCompileShader( handle() );
			GLint st;
			glGetShaderiv( handle(), GL_COMPILE_STATUS, &st );
			if (!st)
			{
				glGetShaderInfoLog( handle(), *len, len, errlog );
				return( false );
			}
			return( true );
		}
		bool info( char *log, GLsizei *len )
		{
			glGetShaderInfoLog( handle(), *len, len, log );
		}
	};
	
	class vertexshader: public shader
	{
	public:
		vertexshader( const char *text ) : shader(text)
		{
			handle( glCreateShader( GL_VERTEX_SHADER ) );
		}
	};

	class geometryshader: public shader
	{
	public:
		geometryshader( const char *text ) : shader(text)
		{
			handle( glCreateShader( GL_GEOMETRY_SHADER_ARB ) );
		}
	};

	class fragmentshader: public shader
	{
	public:
		fragmentshader( const char *text ) : shader(text)
		{
			handle( glCreateShader( GL_FRAGMENT_SHADER ) );
		}
	};

	//
	// interface base classes
	//

	class interface
	{
		const char *_name;
		GLint _binding;
	public:
		interface( const char *name ) : _name(name), _binding(-1) {}
		bool init( GLint index )
		{
			_binding = index;
			if (_binding < 0)
			{
				printf( "glsl: error binding %s\n", _name );
				if (strict) return( false );
			}
			return( true );
		}
		const char *name() const { return( _name ); }
		GLint binding() const { return( _binding ); }
	};

	class interface_in: public interface
	{
	public:
		interface_in( const char *name ) : interface(name) {}
	};

	class interface_uniform: public interface
	{
	public:
		interface_uniform( const char *name ) : interface(name) {}
	};

	#define _GLSL_DECL_1(type) type v0
	#define _GLSL_DECL_2(type) type v0, type v1
	#define _GLSL_DECL_3(type) type v0, type v1, type v2
	#define _GLSL_DECL_4(type) type v0, type v1, type v2, type v3
	#define __GLSL_DECL(name,type) name(type)
	#define _GLSL_DECL(len,type) __GLSL_DECL(_GLSL_DECL_##len,type)
	#define _GLSL_ARG_1 v0
	#define _GLSL_ARG_2 v0, v1
	#define _GLSL_ARG_3 v0, v1, v2
	#define _GLSL_ARG_4 v0, v1, v2, v3
	#define __GLSL_ARG(name) name
	#define _GLSL_ARG(len) __GLSL_ARG(_GLSL_ARG_##len)

	#define _GLSL_LIST_SCALARS(_) \
	_(float,GLfloat,1,f) \
	_(vec2,GLfloat,2,f) \
	_(vec3,GLfloat,3,f) \
	_(vec4,GLfloat,4,f) \
	_(bool,GLint,1,i) \
	_(bvec2,GLint,2,i) \
	_(bvec3,GLint,3,i) \
	_(bvec4,GLint,4,i) \
	_(int,GLint,1,i) \
	_(ivec2,GLint,2,i) \
	_(ivec3,GLint,3,i) \
	_(ivec4,GLint,4,i) \
	_(uint,GLuint,1,ui) \
	_(uvec2,GLuint,2,ui) \
	_(uvec3,GLuint,3,ui) \
	_(uvec4,GLuint,4,ui)

	#define _GLSL_LIST_MATS(_) \
	_(mat2,GLfloat,2,f) \
	_(mat3,GLfloat,3,f) \
	_(mat4,GLfloat,4,f) \
	_(mat2x3,GLfloat,2x3,f) \
	_(mat3x2,GLfloat,3x2,f) \
	_(mat2x4,GLfloat,2x4,f) \
	_(mat4x2,GLfloat,4x2,f) \
	_(mat3x4,GLfloat,3x4,f) \
	_(mat4x3,GLfloat,4x3,f)

	// all this to allow in_* attributes silently drop the normalize for integer types
	#define _VTXATTRIBPTRTYPE_GLfloat(binding,len,ptrtype,normalized,stride,pointer) \
		glVertexAttribPointer(binding,len,ptrtype,normalized,stride,pointer)
	#define _VTXATTRIBPTRTYPE_GLint(binding,len,ptrtype,normalized,stride,pointer) \
		glVertexAttribIPointer(binding,len,ptrtype,stride,pointer)
	#define _VTXATTRIBPTRTYPE_GLuint(binding,len,ptrtype,normalized,stride,pointer) \
		glVertexAttribIPointer(binding,len,ptrtype,stride,pointer)
	#define __MK_VTXATTRIBPTR(type,binding,len,ptrtype,normalized,stride,pointer) \
		type(binding,len,ptrtype,normalized,stride,pointer)
	#define _MK_VTXATTRIBPTR(type,binding,len,ptrtype,normalized,stride,pointer) \
		__MK_VTXATTRIBPTR(type,binding,len,ptrtype,normalized,stride,pointer)
	#define MK_VTXATTRIBPTR(type,binding,len,ptrtype,normalized,stride,pointer) \
		_MK_VTXATTRIBPTR(_VTXATTRIBPTRTYPE_##type,binding,len,ptrtype,normalized,stride,pointer)

	//
	// scalar attributes
	//

	#define _GLSL_MK_IN_SCALAR(name,type,len,t) \
	class in_##name: public interface_in \
	{ \
	public: \
 		in_##name( const char *_name ) : interface_in(_name) {} \
		void bind( GLenum ptrtype, GLsizei stride, const GLvoid *pointer, GLboolean normalized = GL_FALSE ) \
		{ \
			if (binding() >= 0) MK_VTXATTRIBPTR(type, binding(), len, ptrtype, normalized, stride, pointer ); \
		} \
		void enable() { if (binding() >= 0) glEnableVertexAttribArray( binding() ); } \
		void disable() { if (binding() >= 0) glDisableVertexAttribArray( binding() ); } \
	};

	#define _GLSL_MK_IN_SCALAR_ARRAY(name,type,len,t) \
	class in_##name##_array: public interface_in \
	{ \
	public: \
		in_##name##_array( const char *_name ) : interface_in(_name) {} \
		in_##name##_array() : interface_in(NULL) {} \
		void *operator new( size_t s, in_##name##_array *p ) { return p; } \
		void bind( GLenum ptrtype, GLsizei stride, const GLvoid *pointer, GLboolean normalized = GL_FALSE ) \
		{ \
			MK_VTXATTRIBPTR(type, binding(), len, ptrtype, normalized, stride, pointer ); \
		} \
		void enable() { glEnableVertexAttribArray( binding() ); } \
		void disable() { glDisableVertexAttribArray( binding() ); } \
	};

	//
	// matrix attributes
	//

	#define _GLSL_MK_IN_MAT(name,type,size,t) \
	class in_##name: public interface_in \
	{ \
	public: \
		in_##name( const char *_name ) : interface_in(_name) {} \
	};

	#define _GLSL_MK_IN_MAT_ARRAY(name,type,size,t) \
	class in_##name##_array: public interface_in \
	{ \
	public: \
		in_##name##_array( const char *_name ) : interface_in(_name) {} \
		in_##name##_array() : interface_in(NULL) {} \
		void *operator new( size_t s, in_##name##_array *p ) { return p; } \
	};

	//
	// scalar uniforms
	//

	#define _GLSL_MK_UNIFORM_SCALAR(name,type,len,t) \
	class uniform_##name: public interface_uniform \
	{ \
	public: \
		uniform_##name( const char *_name ) : interface_uniform(_name) {} \
		void set( _GLSL_DECL(len,type) ) { glUniform##len##t( binding(), _GLSL_ARG(len) ); } \
		void set( const type *p ) { glUniform##len##t##v( binding(), 1, p ); } \
	};

	#define _GLSL_MK_UNIFORM_SCALAR_ARRAY(name,type,len,t) \
	class uniform_##name##_array: public interface_uniform \
	{ \
	public: \
		uniform_##name##_array( const char *_name ) : interface_uniform(_name) {} \
		uniform_##name##_array() : interface_uniform(NULL) {} \
		void *operator new( size_t s, uniform_##name##_array *p ) { return p; } \
		void set( const type *p, GLsizei count = 1 ) { glUniform##len##t##v( binding(), count, p ); } \
	};

	//
	// matrix uniforms
	//

	#define _GLSL_MK_UNIFORM_MAT(name,type,size,t) \
	class uniform_##name: public interface_uniform \
	{ \
	public: \
		uniform_##name( const char *_name ) : interface_uniform(_name) {} \
		void set( const type *p, GLboolean transpose = GL_FALSE ) { glUniformMatrix##size##t##v( binding(), 1, transpose, p ); } \
	};

	#define _GLSL_MK_UNIFORM_MAT_ARRAY(name,type,size,t) \
	class uniform_##name##_array: public interface_uniform \
	{ \
	public: \
		uniform_##name##_array( const char *_name ) : interface_uniform(_name) {} \
		uniform_##name##_array() : interface_uniform(NULL) {} \
		void *operator new( size_t s, uniform_##name##_array *p ) { return p; } \
		void set( const type *p, GLsizei count = 1, GLboolean transpose = GL_FALSE ) { glUniformMatrix##size##t##v( binding(), count, transpose, p ); } \
	};

	_GLSL_LIST_SCALARS(_GLSL_MK_UNIFORM_SCALAR)
	_GLSL_LIST_SCALARS(_GLSL_MK_UNIFORM_SCALAR_ARRAY)
	_GLSL_LIST_MATS(_GLSL_MK_UNIFORM_MAT)
	_GLSL_LIST_MATS(_GLSL_MK_UNIFORM_MAT_ARRAY)

	_GLSL_LIST_SCALARS(_GLSL_MK_IN_SCALAR)
	_GLSL_LIST_SCALARS(_GLSL_MK_IN_SCALAR_ARRAY)
	_GLSL_LIST_MATS(_GLSL_MK_IN_MAT)
	_GLSL_LIST_MATS(_GLSL_MK_IN_MAT_ARRAY)

	class uniform_sampler: public interface_uniform
	{
	public:
		uniform_sampler( const char *name ) : interface_uniform(name) {}
		void set( GLint index ) { glUniform1i( binding(), index ); }
		void set( const GLint *p ) { glUniform1iv( binding(), 1, p ); }
	};

	class uniform_sampler_array: public interface_uniform
	{
	public:
		uniform_sampler_array( const char *name ) : interface_uniform(name) {}
		void set( const GLint *p, GLsizei count = 1 ) { glUniform1iv( binding(), count, p ); }
	};

} // namespace glsl

#endif // !defined(__glslglue_h)
