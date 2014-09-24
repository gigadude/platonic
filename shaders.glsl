//
// shaders.glsl - sample unified glsl
//

#version 150

//
// various shader samples, generate GLSL:Shader subclasses
//

vertexshader myvtx
{
	in vec3 objvtx[2];
	in vec2 objtex;

	out vec3 gvtx;
	out vec2 gtex;

	uniform float stellate;

	void main( void )
	{
		gvtx = mix( objvtx[0], objvtx[1], stellate );
		gtex = objtex;
	}
}

geometryshader mygeom
{
	layout(triangles) in;
	layout(triangle_strip, max_vertices = 10) out;

	in vec3 gvtx[3];
	in vec2 gtex[3];

	out vec3 opos;
	out vec3 wpos;
	out vec2 wtex;
	flat out vec3 wnrm;
	out vec3 bary;

	// model*view
	uniform mat4 mv;
	// projection
	uniform mat4 p;

	struct
	{
		vec3 opos;
		vec4 clip;
		vec3 wpos;
		vec2 wtex;
		vec3 bary;
	} v[7];
	
	void mk_vtx( int i, vec3 w, float len )
	{
		v[i].opos = gvtx[0] * w.x + gvtx[1] * w.y + gvtx[2] * w.z;
		vec4 m = mv * vec4(v[i].opos, 1);
		v[i].clip = p * m;
		v[i].wpos = m.xyz / m.w;
		v[i].wtex = gtex[0] * w.x + gtex[1] * w.y + gtex[2] * w.z;
		v[i].bary = w;
	}

	void do_vtx( int i )
	{
		gl_Position = v[i].clip;
		opos = v[i].opos;
		wpos = v[i].wpos;
		wtex = v[i].wtex;
		bary = v[i].bary;
		EmitVertex();
	}

	void draw_tri( ivec3 i, int seq )
	{
		vec3 v10 = v[i.y].wpos - v[i.x].wpos;
		vec3 v20 = v[i.z].wpos - v[i.x].wpos;
		wnrm = normalize( cross( -v10, v20 ) );
		if ((seq & 1) == 0) wnrm = -wnrm;
		if (seq == 0)
		{
			do_vtx( i.x );
			do_vtx( i.y );
		}
		do_vtx( i.z );
	}

	void main( void )
	{
		mk_vtx( 0, vec3(1, 0, 0), 0 );
		mk_vtx( 1, vec3(0, 1, 0), 0 );
		mk_vtx( 2, vec3(0, 0, 1), 0 );

		draw_tri( ivec3(0,1,2), 0 );

		EndPrimitive();
	}
}

fragmentshader myfrag
{
	in vec3 opos;
	in vec3 wpos;
	in vec2 wtex;
	flat in vec3 wnrm;
	in vec3 bary;

	uniform vec3 eyepos;
	uniform vec3 lightpos;
	uniform vec3 lightcolor;
	uniform vec3 objcolor;
	uniform vec3 warmcolor;
	uniform vec3 coolcolor;
	uniform float ambient;
	uniform bool outline;
	uniform bool luma;

	const float shininess = 32;
	const float glow = 1.414214;

	out vec4 gl_FragColor;

	void main( void )
	{
		vec3 eyedir = normalize( eyepos - wpos );
		vec3 lightdir = normalize( lightpos - wpos );
		vec3 r = reflect( wnrm, lightdir );
		float l = dot( wnrm, lightdir );
		float k = 0.5 * (l + 1.0);
		vec3 lambert = objcolor * lightcolor * max( ambient, l );
		vec3 gootch = mix( coolcolor, warmcolor, k );
		vec3 spec = lightcolor * pow( max( 0.0, dot( r, eyedir ) ), shininess );
		float b = 1.0;
		if (outline)
		{
			vec3 bcut = 2 * abs( 0.5 - bary );
			b -= smoothstep( 0.95, 1.0, max( max( bcut.x, bcut.y ), bcut.z ) );
		}
		float s = length( opos ) * glow;
		vec3 c = s * s * b * 0.5 * (lambert + gootch) + 0.3 * spec;
		if (luma) c = vec3(c.r * 0.299 + c.g * 0.587 + c.b * 0.114);
		gl_FragColor = vec4( clamp( c, 0, 1 ), 1.0 );
	}
}

//
// myprog - sample shader program linkage; generates a GLSL::Program subclass
//

program myprog
{
	vertexshader myvtx;
	geometryshader mygeom;
	fragmentshader myfrag;
}
