//
// platonic.cpp - show platonic solids evolving into one another
//

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#define GL_GLEXT_PROTOTYPES
#include <GL/freeglut.h>

#include "ssevector.h"

#include "shaders.h"

//
// datatypes
//

// list of polyhedra we can draw/mutate
#define LIST_POLYHEDRA(_) \
_(Tetrahedron) \
_(Cube) \
_(Octahedron) \
_(Dodecahedron) \
_(Icosahedron) \
_(StellaOctangula) \
_(CubeOctahedronCompound) \
_(RhombicDodecahedron) \
_(Cuboctahedron) \
_(DeltoidalIcositetrahedron) \
_(DodecahedronIcosahedronCompound) \
_(RhombicTriacontahedron) \
_(Icosidodecahedron) \
_(DeltoidalHexecontahedron) \

//#define VERBOSE_SIMPLIFY
//#define VERBOSE_TENT

#define WARN(x) do { printf x; } while (0)
#define FAIL(x) do { printf x; exit( 1 ); } while (0)

#define SGN(a,b) ((b > a) - (a > b))

#define CHECK_GLERROR() \
do { \
    GLenum err = glGetError(); \
    if (err != GL_NO_ERROR) \
    	printf( "%s(%d): GL Error %s\n", __FILE__, __LINE__, gluErrorString( err ) ); \
} while (0)

#define MK_SPECIALKEY(x) (0x100|x)
#define SPECIALKEY(x) MK_SPECIALKEY(GLUT_KEY_##x)

// list of menu label, value, case statement and command to invoke
#define LIST_COMMANDS(_) \
_("Animate  [a]",'a',case 'a':,(bAnimate ^= true)) \
_("Outline  [o]",'o',case 'o':,(bOutline ^= true)) \
_("Switch  [n]",'n',case 'n':,(bSwitchPolys ^= true)) \
_("Random  [r]",'r',case 'r':,(bRandom ^= true)) \
_("Red/Cyan  [d]",'d',case 'd':,(bStereo ^= true)) \
_("Tetrahedron  [t]",'t',case 't':,setPolyMode( Tetrahedron )) \
_("Cube  [c]",'c',case 'c':,setPolyMode( Cube )) \
_("Add Pyramids  [p]",'p',case 'p':,setChangeMode( ADD_PYRAMIDS )) \
_("Add Stellations  [s]",'s',case 's':,setChangeMode( ADD_STELLATIONS )) \
_("Simplify  [S]",'S',case 'S':,toggleSimplify()) \
_("Keep All  [0]",'0',case '0':,setKeepMode( KEEP_ALL )) \
_("Keep Points  [1]",'1',case '1':,setKeepMode( KEEP_POINTS )) \
_("Keep Saddles  [2]",'2',case '2':,setKeepMode( KEEP_SADDLES )) \
_("Keep Stellations  [3]",'3',case '3':,setKeepMode( KEEP_STELLATIONS )) \
_("Current Poly  [ ]",' ',case ' ':,selectPoly( 0 )) \
_("Next Poly  []]",']',case ']':,selectPoly( 1 )) \
_("Prev Poly  [[]",'[',case '[':,selectPoly( -1 )) \
_("Exit  [Esc]",27,case 27:,exit(0))

const float M_PHI = 0.5f * (sqrt( 5.0f ) - 1.0f);
const float EPSILON = 1e-6f;

//
// vtx - describes a vertex in the polyhedron
//

class vtx
{
public:
    float	p[2][3];
    float	t[2];
    enum vtype { POINT, STELLATE, SADDLE } type;

    void zero()
    {
	for (int i = 0; i < 2; ++i)
	    for (int j = 0; j < 3; ++j)
		p[i][j] = 0;
	for (int i = 0; i < 2; ++i)
	    t[i] = 0;
	type = POINT;
    }
    void swap()
    {
	for (int j = 0; j < 3; ++j)
	{
	    float tmp = p[0][j];
	    p[0][j] = p[1][j];
	    p[1][j] = tmp;
	}
    }
    vtx &operator =( const vtx &v )
    {
	for (int i = 0; i < 2; ++i)
	    for (int j = 0; j < 3; ++j)
		p[i][j] = v.p[i][j];
	for (int i = 0; i < 2; ++i)
	    t[i] = v.t[i];
	type = v.type;
	return( *this );
    }
    vtx &operator +=( const vtx &v )
    {
	for (int i = 0; i < 2; ++i)
	    for (int j = 0; j < 3; ++j)
		p[i][j] += v.p[i][j];
	for (int i = 0; i < 2; ++i)
	    t[i] += v.t[i];
	return( *this );
    }
    vtx &operator -=( const vtx &v )
    {
	for (int i = 0; i < 2; ++i)
	    for (int j = 0; j < 3; ++j)
		p[i][j] -= v.p[i][j];
	for (int i = 0; i < 2; ++i)
	    t[i] -= v.t[i];
	return( *this );
    }
    vtx &operator *=( float s )
    {
	for (int i = 0; i < 2; ++i)
	    for (int j = 0; j < 3; ++j)
		p[i][j] *= s;
	for (int i = 0; i < 2; ++i)
	    t[i] *= s;
	return( *this );
    }
    vtx &operator *=( const float s[] )
    {
	for (int i = 0; i < 2; ++i)
	    for (int j = 0; j < 3; ++j)
		p[i][j] *= s[i];
	for (int i = 0; i < 2; ++i)
	    t[i] *= s[i];
	return( *this );
    }
    float dot( const vtx &v, int fromvtx = 0 ) const
    {
	float d =
	    p[fromvtx][0] * v.p[fromvtx][0] + 
	    p[fromvtx][1] * v.p[fromvtx][1] +
	    p[fromvtx][2] * v.p[fromvtx][2];
	return( d );
    }
    float length( int fromvtx = 0 ) const
    {
	return( sqrtf( dot( *this, fromvtx ) ) );
    }
    float normalize( int fromvtx = 0 )
    {
	float len = length( fromvtx );
	if (len > 0) *this *= 1.0f / len;
	return( len );
    }
    float distance( const vtx &v, int fromvtx = 0 ) const
    {
	return( sqrtf( dot( v, fromvtx ) ) );
    }
    vtx cross( const vtx &v ) const
    {
	vtx t;
	for (int i = 0; i < 2; ++i)
	{
	    t.p[i][0] = p[i][1] * v.p[i][2] - v.p[i][1] * p[i][2];
	    t.p[i][1] = p[i][2] * v.p[i][0] - v.p[i][2] * p[i][0];
	    t.p[i][2] = p[i][0] * v.p[i][1] - v.p[i][0] * p[i][1];
	}
	return( t );
    }
    void setType( vtype t ) { type = t; }
    bool isPoint() const { return( type == POINT ); }
    bool isStellate() const { return( type == STELLATE ); }
    bool isSaddle() const { return( type == SADDLE ); }
};

// offsets to fields once the vbo is bound (offsetof is broken for gcc, doesn't handle arrays w/ constant indexes)
#define vbof(t,f) (size_t(&(static_cast<t *>((void *)1)->f)) - 1)
#define VBO_VTX_OFFSET(fld)		((GLvoid *)vbof(vtx,fld))
#define VBO_VTX_OFFSET_ARRAY(fld,idx)	((GLvoid *)(vbof(vtx,fld[0]) + (idx)*(vbof(vtx,fld[1]) - vbof(vtx,fld[0]))))

//
// face - describes a face in a polyhedron
//

static const int MAX_FACE_VERTS = 15;

class face
{
public:
    int		verts;			// number of vertexes in face
    int		vi[MAX_FACE_VERTS];	// indexes of vtx structs
    face() { verts = 0; }
    face( const face &copy ) :
	verts(copy.verts)
    {
	for (int i = 0; i < copy.verts; ++i) vi[i] = copy.vi[i];
    }
};

//
// edgemap - map directed [fromvert][tovert] edges to an index
//

class edgemap
{
    int		verts;
    int		*e;			// [fromvert][tovert] matrix mapping edge -> int
    int edgeidx( int fromvi, int tovi ) const
    {
	return( fromvi * verts + tovi );
    }
public:
    edgemap( int verts ) :
	verts(verts)
    {
	e = new int[verts * verts];
	clear();
    }
    edgemap( const edgemap &copy ) :
	verts(copy.verts)
    {
	e = new int[verts * verts];
	for (int i = 0; i < verts * verts; ++i) e[i] = copy.e[i];
    }
    ~edgemap()
    {
	delete [] e;
	verts = 0;
    }
    void clear()
    {
	for (int i = 0; i < verts * verts; ++i) e[i] = -1;
    }
    void clear( int fromvi, int tovi )
    {
	e[edgeidx( fromvi, tovi )] = -1;
    }
    int get( int fromvi, int tovi ) const
    {
	return( e[edgeidx( fromvi, tovi )] );
    }
    void set( int fromvi, int tovi, int val )
    {
	e[edgeidx( fromvi, tovi )] = val;
    }
};

//
// poly - polyhedron descriptor
//

class poly
{
public:
    int		verts, faces;
    int		next_face;
    edgemap	e;			// [fromvert][tovert] matrix mapping edge -> faceidx
    vtx		*v;
    face	*f;

    poly( int verts, int faces ) :
	verts(verts),
	faces(faces),
	next_face(0),
	e(verts)
    {
	e.clear();
	v = new vtx[verts];
	f = new face[faces];
    }
    poly( const poly &copy ) :
	verts(copy.verts),
	faces(copy.next_face),
	next_face(0),
	e(copy.verts)
    {
	e.clear();
	v = new vtx[verts];
	f = new face[faces];
	for (int i = 0; i < verts; ++i) setvtx( i, copy.v[i] );
	for (int i = 0; i < copy.next_face; ++i) addface( copy.f[i].verts, copy.f[i].vi );
    }
    ~poly()
    {
	delete[] v;
	delete[] f;
    }
    int setvtx( int vi, const vtx &v0 )
    {
	if (vi >= verts) return( -1 );
	v[vi] = v0;
	return( vi );
    }
    void setedge( int fromvi, int tovi, int fi )
    {
	e.set( fromvi, tovi, fi );
    }
    int edges() const
    {
	// for convex polyhedra, verts - edges + faces = 2
	return( verts + faces - 2 );
    }
    int faceidx( int fromvi, int tovi ) const
    {
	return( e.get( fromvi, tovi ) );
    }
    int addface( int face_verts, const int idx[] )
    {
	if (next_face >= faces) return( -1 );
	int fi = next_face++;
	f[fi].verts = face_verts;
	int last_vi = idx[face_verts - 1];
	for (int i = 0; i < face_verts; ++i)
	{
	    f[fi].vi[i] = idx[i];
	    setedge( last_vi, idx[i], fi );
	    last_vi = idx[i];
	}
	return( fi );
    }
    void clearfaces()
    {
	e.clear();
	next_face = 0;
    }
    void normalize( int i )
    {
	float maxlen = 0;
	for (int j = 0; j < verts; ++j)
	{
	    float len = v[j].dot( v[j], i );
	    if (len > maxlen) maxlen = len;
	}
	float scale = 1.0f / sqrtf( maxlen );
	for (int j = 0; j < verts; ++j)
	{
	    v[j].p[i][0] *= scale;
	    v[j].p[i][1] *= scale;
	    v[j].p[i][2] *= scale;
	}
    }
    vtx facecenter( int i ) const
    {
	vtx t = v[f[i].vi[0]];
	for (int j = 1; j < f[i].verts; ++j)
	    t += v[f[i].vi[j]];
	float s = 1.0f / f[i].verts;
	t *= s;
	return( t );
    }
    vtx edgecenter( int fromvi, int tovi ) const
    {
	vtx t = v[fromvi];
	t += v[tovi];
	t *= 0.5f;
	return( t );
    }
};

//
// scaleToNeighbors - scale a face-centered point to meet the plane containing an edge most perpendicular to the origin
//

float scaleToNeighbors( const poly *pin, int tovtx, int i, const vtx &fc )
{
    const face &iface = pin->f[i];
    const vtx &p0 = pin->v[iface.vi[0]];
    const vtx &p1 = pin->v[iface.vi[1]];
    vtx pd = p1;
    pd -= p0;
    vtx a = p1.cross( p0 );
    vtx n = a.cross( pd );
    n.normalize( tovtx );
    float edgedist = p0.dot( n, tovtx );
    float facedist = fc.dot( n, tovtx );
    return( edgedist / facedist );
}

//
// addPyramids - generate a vertex in each face center and make a pyramid using vertex positions
//

poly *addPyramids( const poly *pin, int tovtx, float scale = 0.0f )
{
    int verts = pin->verts + pin->faces;
    int faces = 0;
    for (int i = 0; i < pin->faces; ++i)
	faces += pin->f[i].verts;

    poly *pout = new poly( verts, faces );

    // copy the original verts
    for (int i = 0; i < pin->verts; ++i)
	pout->setvtx( i, pin->v[i] );

    // add a new vertex in each face center
    for (int i = 0; i < pin->faces; ++i)
    {
	vtx fc = pin->facecenter( i );
	float flen = (scale > 0.0f) ? scale : scaleToNeighbors( pin, tovtx, i, fc );
	fc.p[tovtx][0] *= flen;
	fc.p[tovtx][1] *= flen;
	fc.p[tovtx][2] *= flen;
	fc.setType( vtx::STELLATE );
	pout->setvtx( pin->verts + i, fc );
    }

    // generate a pyramid of faces using each original face and the center vertex
    for (int i = 0; i < pin->faces; ++i)
    {
	int last_vi = pin->f[i].vi[pin->f[i].verts - 1];
	for (int j = 0; j < pin->f[i].verts; ++j)
	{
	    int idx[3];
	    idx[0] = last_vi;
	    idx[1] = pin->f[i].vi[j];
	    idx[2] = pin->verts + i;
	    pout->addface( 3, idx );
	    last_vi = pin->f[i].vi[j];
	}
    }

    pout->normalize( 0 );
    pout->normalize( 1 );
    
    return( pout );
}

//
// addStellations - generate a vertex in each face center and make a pyramid using edge midpoints
//

poly *addStellations( const poly *pin, int tovtx, float scale = 0.0f )
{
    //printf( "addStellations pin verts = %d, faces = %d, edges = %d\n", pin->verts, pin->faces, pin->edges() );
    int verts = pin->verts + pin->faces + pin->edges();
    int faces = 0;
    for (int i = 0; i < pin->faces; ++i)
	faces += 2 * pin->f[i].verts;

    poly *pout = new poly( verts, faces );

    // copy the original verts
    for (int i = 0; i < pin->verts; ++i)
	pout->setvtx( i, pin->v[i] );

    // add a new vertex in each face center
    for (int i = 0; i < pin->faces; ++i)
    {
	vtx fc = pin->facecenter( i );
	float flen = (scale > 0.0f) ? scale : scaleToNeighbors( pin, tovtx, i, fc );
	fc.p[tovtx][0] *= flen;
	fc.p[tovtx][1] *= flen;
	fc.p[tovtx][2] *= flen;
	fc.setType( vtx::STELLATE );
	pout->setvtx( pin->verts + i, fc );
    }

    // add each edge midpoint
    edgemap midvtx( pin->verts );
    int next_mid = pin->verts + pin->faces;
    for (int i = 0; i < pin->faces; ++i)
    {
	int last_vi = pin->f[i].vi[pin->f[i].verts - 1];
	for (int j = 0; j < pin->f[i].verts; ++j)
	{
	    int this_vi = pin->f[i].vi[j];
	    if (midvtx.get( last_vi, this_vi ) < 0)
	    {
		if (next_mid >= verts) FAIL(( "addStellations overflow [%d]!", next_mid ));
		// map both directions so each edge shares one vertex
		midvtx.set( last_vi, this_vi, next_mid );
		midvtx.set( this_vi, last_vi, next_mid );
		//printf( "add midpoint %d (%d,%d)\n", next_mid, last_vi, this_vi );
		vtx t = pin->v[last_vi];
		t += pin->v[this_vi];
		t *= 0.5f;
		t.setType( vtx::SADDLE );
		pout->setvtx( next_mid, t );
		++next_mid;
	    }
	    last_vi = this_vi;
	}
    }
    //printf( "added %d midpoints for %d edges\n", next_mid - (pin->verts + pin->faces), pin->edges() );

    // stellate the face using the midpoints and the center
    for (int i = 0; i < pin->faces; ++i)
    {
	//printf( "adding face %d\n", i );
	int last_vi = pin->f[i].vi[pin->f[i].verts - 1];
	int center_vi = pin->verts + i;
	for (int j = 0; j < pin->f[i].verts; ++j)
	{
	    int this_vi = pin->f[i].vi[j];
	    int next_vi = pin->f[i].vi[((j + 1) < pin->f[i].verts) ? (j + 1) : 0];
	    int last_mid_vi = midvtx.get( last_vi, this_vi );
	    int next_mid_vi = midvtx.get( this_vi, next_vi );

	    int idx[3];

	    idx[0] = last_mid_vi;
	    idx[1] = next_mid_vi;
	    idx[2] = center_vi;
	    pout->addface( 3, idx );
	    //printf( "(%d,%d,%d)\n", last_mid_vi, next_mid_vi, center_vi );

	    idx[0] = last_mid_vi;
	    idx[1] = this_vi;
	    idx[2] = next_mid_vi;
	    pout->addface( 3, idx );
	    //printf( "(%d,%d,%d)\n", last_mid_vi, this_vi, next_mid_vi );

	    last_vi = this_vi;
	}
    }

    pout->normalize( 0 );
    pout->normalize( 1 );
    
    return( pout );
}

//
// addTents - add tents to the faces of a simplified cube transforming it to both a dodecahedron and a rhombic dodecahedron
//

void colorCubeTets( const poly *pin, int i, int *vtet )
{
    for (int j = 0 + (i == 0); j < pin->verts; j += 1 + (j == i))
    {
	if ((vtet[j] == -1) && (pin->e.get( i, j ) >= 0))
	{
	    vtet[j] = vtet[i] ^ 1;
	    colorCubeTets( pin, j, vtet );
	}
    }
}

poly *addTents( const poly *pin, int tovtx, float fromscale, float toscale )
{
    int fromvtx = tovtx ^ 1;

    // color the cube verts based on which tetrahedron they belong to
    int *vtet = new int[pin->verts];
    for (int i = 0; i < pin->verts; ++i) vtet[i] = -1;
    vtet[0] = 0;
    colorCubeTets( pin, 0, vtet );

    #if defined(VERBOSE_TENT)
	printf( "addTents: pin( %d, %d ), %d->%d (%f->%f), vtet =", pin->verts, pin->next_face, fromvtx, tovtx, fromscale, toscale );
	for (int i = 0; i < pin->verts; ++i) printf( " %2d", vtet[i] );
	printf( "\n" );
    #endif

    poly *pout = new poly( pin->verts + pin->next_face * 2, pin->next_face * 4 );

    // copy the original verts
    for (int i = 0; i < pin->verts; ++i)
	pout->setvtx( i, pin->v[i] );

    // scalars, shp1 = 1 + h, s1mh2 = 1 - h^2 where h is the from/to scale
    float shp1[2];
    shp1[fromvtx] = 1.0f + fromscale;
    shp1[tovtx] = 1.0f + toscale;
    float s1mh2[2];
    s1mh2[fromvtx] = 1.0f - fromscale * fromscale;
    s1mh2[tovtx] = 1.0f - toscale * toscale;

    // generate the tent verts and faces for each original face
    for (int i = 0; i < pin->next_face; ++i)
    {
	face &iface = pin->f[i];

	#if defined(VERBOSE_TENT)
	    printf( "addTent: face %d ([%d]", i, iface.verts );
	    for (int j = 0; j < iface.verts; ++j) printf( " %d/%d", iface.vi[j], vtet[iface.vi[j]] );
	    printf( ")\n" );
	#endif

	int j = 0;
	int last_vi = iface.verts - 1;

	// rotate face verts if j==0 is for the other tetrahedron in the cube
	if (vtet[iface.vi[0]])
	{
	    j = 1;
	    last_vi = 0;
	}

	// faces look like:
	//       t0
	//  lvi--+-->j
	//   ^   |   |
	//   |   |   V
	//  j+2<-+--j+1
	//       t1

	vtx fc = pin->facecenter( i );
	vtx t0 = pin->edgecenter( iface.vi[last_vi], iface.vi[j] );
	vtx t1 = pin->edgecenter( iface.vi[j + 1], iface.vi[j + 2] );

	// top center of the tent is face center scaled by 1 + h
	vtx tc = fc;
	tc *= shp1;

	// directions along tent top are distance to midpoint scaled by 1 - h^2
	t0 -= fc; t0 *= s1mh2; t0 += tc;
	t1 -= fc; t1 *= s1mh2; t1 += tc;

	// store the verts indexed by face # * 2 after the original verts
	int tvi0 = pin->verts + i * 2;
	int tvi1 = tvi0 + 1;
	#if defined(VERBOSE_TENT)
	    printf( "adding vtx %d/%d\n", tvi0, tvi1 );
	#endif
	pout->setvtx( tvi0, t0 );
	pout->setvtx( tvi1, t1 );

	// build the new faces
	int idx[4];

	idx[0] = iface.vi[j];
	idx[1] = tvi0;
	idx[2] = iface.vi[last_vi];
	pout->addface( 3, idx );

	idx[0] = iface.vi[j];
	idx[1] = iface.vi[j + 1];
	idx[2] = tvi1;
	idx[3] = tvi0;
	pout->addface( 4, idx );
	
	idx[0] = iface.vi[j + 1];
	idx[1] = iface.vi[j + 2];
	idx[2] = tvi1;
	pout->addface( 3, idx );

	idx[0] = iface.vi[j + 2];
	idx[1] = iface.vi[last_vi];
	idx[2] = tvi0;
	idx[3] = tvi1;
	pout->addface( 4, idx );
    }

    delete [] vtet;

    #if defined(VERBOSE_TENT)
	printf( "addTents: pout %d verts %d faces\n", pout->verts, pout->next_face );
    #endif

    pout->normalize( 0 );
    pout->normalize( 1 );

    return( pout );
}

//
// removeTents - remove tents from the verts of a simplified octahedron transforming it to a icosahedron
//

void colorTentFace( poly *pout, const poly *pin, int tovtx, int fi, int parity, edgemap &phivtx, int &nextphivtx )
{
    const face &iface = pin->f[fi];
    int last_vi = iface.verts - 1;
    for (int i = 0; i < iface.verts; ++i)
    {
	int vi0 = iface.vi[last_vi];
	int vi1 = iface.vi[i];
	if (phivtx.get( vi0, vi1 ) < 0)
	{
	    vtx v0 = pin->v[vi0];
	    vtx v1 = pin->v[vi1];
	    if (parity) v1 *= M_PHI;
	    else v0 *= M_PHI;
	    vtx v = v0;
	    v += v1;
	    v *= 1.0f / (1.0f + M_PHI);
	    pout->setvtx( nextphivtx, v );
	    phivtx.set( vi0, vi1, nextphivtx );
	    phivtx.set( vi1, vi0, nextphivtx );
	    ++nextphivtx;
	    colorTentFace( pout, pin, tovtx, pin->e.get( vi1, vi0 ), parity ^ 1, phivtx, nextphivtx );
	}
	last_vi = i;
    }
}

poly *removeTents( const poly *pin, int tovtx )
{
    poly *pout = new poly( pin->verts + pin->edges(), pin->next_face * 4 );

    // generate the bottom of the tents
    int fromvtx = tovtx ^ 1;
    for (int i = 0; i < pin->verts; ++i)
    {
	vtx v = pin->v[i];
	v.p[tovtx][0] = v.p[fromvtx][0] * M_PHI;
	v.p[tovtx][1] = v.p[fromvtx][1] * M_PHI;
	v.p[tovtx][2] = v.p[fromvtx][2] * M_PHI;
	pout->setvtx( i, v );
    }

    // generate the midpoints by subdividing edges by phi
    int nextphivtx = pin->verts;
    edgemap phivtx( pin->verts );
    colorTentFace( pout, pin, tovtx, 0, 1, phivtx, nextphivtx );

    for (int i = 0; i < pin->next_face; ++i)
    {
	face &iface = pin->f[i];

	int idx[3];

	idx[0] = phivtx.get( iface.vi[2], iface.vi[0] );
	idx[1] = phivtx.get( iface.vi[0], iface.vi[1] );
	idx[2] = phivtx.get( iface.vi[1], iface.vi[2] );
	pout->addface( 3, idx );

	idx[0] = iface.vi[0];
	idx[1] = phivtx.get( iface.vi[0], iface.vi[1] );
	idx[2] = phivtx.get( iface.vi[2], iface.vi[0] );
	pout->addface( 3, idx );

	idx[0] = iface.vi[1];
	idx[1] = phivtx.get( iface.vi[1], iface.vi[2] );
	idx[2] = phivtx.get( iface.vi[0], iface.vi[1] );
	pout->addface( 3, idx );

	idx[0] = iface.vi[2];
	idx[1] = phivtx.get( iface.vi[2], iface.vi[0] );
	idx[2] = phivtx.get( iface.vi[1], iface.vi[2] );
	pout->addface( 3, idx );
    }

    #if defined(VERBOSE_TENT)
	printf( "removeTents: pout %d verts %d faces\n", pout->verts, pout->next_face );
    #endif

    pout->normalize( 0 );
    pout->normalize( 1 );

    return( pout );
}

//
// keepVertexType - keep the specified verts as-is and make everything else morph to co-planar with those verts
//

poly *keepVertexType( const poly *pin, int tovtx, vtx::vtype type )
{
    poly *pout = new poly( *pin );
    int fromvtx = tovtx ^ 1;

    // check each vertex
    for (int i = 0; i < pin->verts; ++i)
    {
	// copy the start to final vertex
	vtx t = pin->v[i];
	t.p[tovtx][0] = t.p[fromvtx][0];
	t.p[tovtx][1] = t.p[fromvtx][1];
	t.p[tovtx][2] = t.p[fromvtx][2];
	
	if (t.type != type)
	{
	    // saddles are already in the right place if we're keeping POINT/STELLATE verts
	    if (t.type != vtx::SADDLE)
	    {
		int vert_edges = 0;
		for (int j = 0; j < pin->verts; ++j)
		{
		    if (j == i) continue;
		    if (pin->e.get( i, j ) >= 0)
		    {
			if (vert_edges++ == 0)
			{
			    t.p[tovtx][0] = pin->v[j].p[fromvtx][0];
			    t.p[tovtx][1] = pin->v[j].p[fromvtx][1];
			    t.p[tovtx][2] = pin->v[j].p[fromvtx][2];
			}
			else
			{
			    t.p[tovtx][0] += pin->v[j].p[fromvtx][0];
			    t.p[tovtx][1] += pin->v[j].p[fromvtx][1];
			    t.p[tovtx][2] += pin->v[j].p[fromvtx][2];
			}
		    }
		}
		if (vert_edges > 1 )
		{
		    float s = 1.0f / vert_edges;
		    t.p[tovtx][0] *= s;
		    t.p[tovtx][1] *= s;
		    t.p[tovtx][2] *= s;
		}
		t.setType( vtx::STELLATE );
	    }
	}
	else t.setType( vtx::POINT ); // keepers become the new points

	pout->setvtx( i, t );
    }

    pout->normalize( 0 );
    pout->normalize( 1 );
    
    return( pout );
}

//
// simplify - merge co-planar faces
//

void simplifyFace( poly *pout, face &oface,
    const poly *pin, int &last_vtx, int i,
    const edgemap &flat, const int *vmap, int *fmap )
{
     if (fmap[i] >= 0) return; // already merged
    fmap[i] = pout->next_face;
    const face &iface = pin->f[i];

    #if defined(VERBOSE_SIMPLIFY)
	printf( "simplifyFace: %d ([%d]", i, iface.verts );
	for (int j = 0; j < iface.verts; ++j) printf( " %d", iface.vi[j] );
	printf( ") @ %d -> %d ([%d]", last_vtx, pout->next_face, oface.verts );
	for (int j = 0; j < oface.verts; ++j) printf( " %d", oface.vi[j] );
	printf( ")\n" );
    #endif

    // find starting face vertex
    int j = 0;
    int last_vi = iface.verts - 1;
    if ((last_vtx >= 0) && (iface.vi[last_vi] != last_vtx))
    {
	// we checked j==0 above, should find it at 1..verts-1
	for (j = 1; j < iface.verts; ++j)
	{
	    last_vi = j - 1;
	    if (iface.vi[last_vi] == last_vtx)
		break;
	}
	if (j >= iface.verts) FAIL(( "simplifyFace: face %d vtx %d not found!\n", i, last_vtx ));
    }
    else last_vtx = iface.vi[last_vi];

    for (int k = 0; k < iface.verts; ++k)
    {
	int fromvi = iface.vi[last_vi];
	int tovi = iface.vi[j];

	#if defined(VERBOSE_SIMPLIFY)
	    printf( "checking edge %d-%d\n", fromvi, tovi );
	#endif

	// if this edge is flat, merge the faces
	if (flat.get( fromvi, tovi ) == 1)
	{
	    // get the index of the other face and merge it starting at last_vtx
	    #if defined(VERBOSE_SIMPLIFY)
		printf( "checking opposite face %d (this is %d)\n", pin->faceidx( tovi, fromvi ), pin->faceidx( fromvi, tovi ) );
	    #endif
	    simplifyFace( pout, oface, pin, last_vtx, pin->faceidx( tovi, fromvi ), flat, vmap, fmap );
	}
	else if (vmap[tovi] >= 0)
	{
	    // add the vertex if we're keeping it
	    #if defined(VERBOSE_SIMPLIFY)
		printf( "adding vtx %d -> %d\n", tovi, vmap[tovi] );
	    #endif
	    oface.vi[oface.verts++] = vmap[tovi];
	}

	last_vtx = tovi;
	if (++j >= iface.verts) j = 0;
	if (++last_vi >= iface.verts) last_vi = 0;
    }
    #if defined(VERBOSE_SIMPLIFY)
	printf( "simplifyFace finished %d (last_vtx = %d)\n", i, last_vtx );
    #endif
}

#define PVTX(v) printf( "%s (%f,%f,%f)\n", #v, v.p[fromvtx][0], v.p[fromvtx][1], v.p[fromvtx][2] )

poly *simplify( const poly *pin, int fromvtx, bool convex = false )
{
    #if defined(VERBOSE_SIMPLIFY)
	printf( "simplify %d verts, %d faces [%d]\n", pin->verts, pin->next_face, fromvtx );
    #endif

    // find all flat edges
    edgemap flat( pin->verts );
    int ekept = 0;
    for (int i = 0; i < pin->verts; ++i)
    {
	for (int j = i + 1; j < pin->verts; ++j)
	{
	    int f0 = pin->faceidx( i, j );
	    int f1 = pin->faceidx( j, i );
	    if ((f0 < 0) || (f1 < 0)) continue;

	    vtx ec = pin->edgecenter( i, j );
	    vtx fc0 = pin->facecenter( f0 );
	    vtx fc1 = pin->facecenter( f1 );
	    vtx f00 = pin->v[pin->f[f0].vi[0]];
	    vtx f01 = pin->v[pin->f[f0].vi[1]];
	    vtx f10 = pin->v[pin->f[f1].vi[0]];
	    vtx f11 = pin->v[pin->f[f1].vi[1]];
	    f00 -= fc0;
	    f01 -= fc0;
	    f10 -= fc1;
	    f11 -= fc1;
	    vtx n0 = f00.cross( f01 );
	    vtx n1 = f10.cross( f11 );
	    // the length of co-linear normals dotted is just the lengths multiplied
	    float n0len = n0.length( fromvtx );
	    float n1len = n1.length( fromvtx );
	    float ndotlen = n0len * n1len;
	    float diff = (ndotlen > 0) ? fabsf( ndotlen - n0.dot( n1, fromvtx ) ) : 1.0f;

	    #if defined(VERBOSE_SIMPLIFY)
		printf( "diff edge %d-%d (face %d/%d) = %f ", i, j, f0, f1, diff );
	    #endif
	    if (diff < EPSILON)
	    {
		flat.set( i, j, 1 );
		flat.set( j, i, 1 );
		#if defined(VERBOSE_SIMPLIFY)
		    printf( "flat\n" );
		#endif
	    }
	    else
	    {
		flat.set( i, j, 0 );
		flat.set( j, i, 0 );
		++ekept;
		#if defined(VERBOSE_SIMPLIFY)
		    printf( "keep %d\n", ekept );
		#endif
	    }
	}
    }

    // eliminate all verts which only connect to flat edges
    int *vmap = new int[pin->verts];
    int *vkeep = new int[pin->verts];
    int vkept = 0;
    for (int i = 0; i < pin->verts; ++i)
    {
	vmap[i] = -1;
	// count nonflat edges to detect colinear verts
	int dihedralcount = 0;
	int evi[2];

	for (int j = 0 + (i == 0); j < pin->verts; j += 1 + (j == i))
	{
	    if (flat.get( i, j ) == 0) // valid and not flat
	    {
		#if defined(VERBOSE_SIMPLIFY)
		    printf( "edge %d-%d not flat\n", i, j );
		#endif
		if (dihedralcount < 2) evi[dihedralcount] = j;
		++dihedralcount;
	    }
	}
	if (dihedralcount > 0)
	{
	    bool keep = true;
	    // check if we've already kept another copy of this vert
	    for (int j = 0; j < i; ++j)
	    {
		vtx dij = pin->v[j];
		dij -= pin->v[i];
		float d = dij.dot( dij, fromvtx );
		if (d < EPSILON)
		{
		    #if defined(VERBOSE_SIMPLIFY)
			printf( "vtx %d duplicates %d\n", i, j );
		    #endif
		    vmap[i] = vmap[j];
		    keep = false;
		    break;
		}
	    }
	    if (keep && (dihedralcount == 2))
	    {
		// check for colinearality
		vtx v0 = pin->v[evi[0]];
		vtx v1 = pin->v[evi[1]];
		v0 -= pin->v[i];
		v0.normalize();
		v1 -= pin->v[i];
		v1.normalize();
		vtx d = v0.cross( v1 );
		if (d.dot( d, fromvtx ) < EPSILON)
		{
		    keep = false;
		    // erase the two contributing edges and mark the new edge as non-flat
		    flat.set( i, evi[0], -1 );
		    flat.set( evi[0], i, -1 );
		    flat.set( i, evi[1], -1 );
		    flat.set( evi[1], i, -1 );
		    flat.set( evi[0], evi[1], 0 );
		    flat.set( evi[1], evi[0], 0 );
		    --ekept;
		    #if defined(VERBOSE_SIMPLIFY)
			printf( "vtx %d colinear with %d-%d (%f)\n", i, evi[0], evi[1], d.dot( d, fromvtx ) );
			PVTX(pin->v[evi[0]]);
			PVTX(pin->v[i]);
			PVTX(pin->v[evi[1]]);
			PVTX(v0);
			PVTX(v1);
			PVTX(d);
		    #endif
		}
	    }
	    if (keep)
	    {
		#if defined(VERBOSE_SIMPLIFY)
		    printf( "vtx %d -> %d\n", i, vkept );
		#endif
		vmap[i] = vkept;
		vkeep[vkept++] = i;
	    }
	}
    }

    // we now know enough to allocate the new polyhedron
    // for convex polyhedra, verts - edges + faces = 2
    // faces = 2 - verts + edges
    #if defined(VERBOSE_SIMPLIFY)
	printf( "simplify: kept %d edges, %d verts from (%d,%d)\n\n", ekept, vkept, pin->verts, pin->edges() );
    #endif
    poly *pout = new poly( vkept, 2 - vkept + ekept );
    int tovtx = fromvtx ^ 1;
    for (int i = 0; i < vkept; ++i)
    {
	// copy the start to final vertex
	vtx t = pin->v[vkeep[i]];
	t.p[tovtx][0] = t.p[fromvtx][0];
	t.p[tovtx][1] = t.p[fromvtx][1];
	t.p[tovtx][2] = t.p[fromvtx][2];
	if (convex) t.type = vtx::POINT;
	pout->setvtx( i, t );
    }

    // check edges on each face, merge faces which are coplanar
    int *fmap = new int[pin->next_face];
    for (int i = 0; i < pin->next_face; ++i) fmap[i] = -1;

    for (int i = 0; i < pin->next_face; ++i)
    {
	face f;
	int last_vtx = -1;
	simplifyFace( pout, f, pin, last_vtx, i, flat, vmap, fmap );
	if (f.verts > 0)
	{
	    #if defined(VERBOSE_SIMPLIFY)
		printf( "adding face %d ([%d]", pout->next_face, f.verts );
		for (int j = 0; j < f.verts; ++j) printf( " %d", f.vi[j] );
		printf( ")\n\n" );
	    #endif
	    pout->addface( f.verts, f.vi );
	}
    }

    pout->normalize( 0 );
    pout->normalize( 1 );

    #if defined(VERBOSE_SIMPLIFY)
	printf( "\nsimplify: %d faces\n", pout->next_face );
    #endif

    delete [] vmap;
    delete [] vkeep;
    delete [] fmap;
    
    return( pout );
}

//
// flip - create a new poly with the opposite morph direction
//

poly *flip( const poly *pin )
{
    poly *pout = new poly( *pin );
    for (int i = 0; i < pout->verts; ++i) pout->v[i].swap();
    return( pout );
}

//
// drawobj - manages vertex buffer objects for a drawable object
//

class drawobj
{
    enum	{ VTX, IDX, NUM_VBOS };
    GLuint 	vbo[2];		// 0 == verts, 1 == indexes

public:
    int		verts;
    int		tris;

    drawobj() : verts(0), tris(0)
    {
	glGenBuffers( NUM_VBOS, vbo );
	if (!vbo[VTX] || !vbo[IDX]) FAIL(( "Failed to create VBOs!" ));
    }
    ~drawobj()
    {
	glDeleteBuffers( NUM_VBOS, vbo );
    }
    void set( poly *p )
    {
	verts = p->verts;
	glBindBuffer( GL_ARRAY_BUFFER, vbo[VTX] );
	glBufferData( GL_ARRAY_BUFFER, sizeof(p->v[0]) * p->verts, &p->v[0], GL_STATIC_DRAW );
	//glBufferSubData( GL_ARRAY_BUFFER, 0, sizeof(v), &v[0] );
	glBindBuffer( GL_ARRAY_BUFFER, 0 );

	tris = 0;
	for (int fi = 0; fi < p->faces; ++fi)
	    tris += p->f[fi].verts - 2;

	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, vbo[IDX] );
	glBufferData( GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * 3 * tris, NULL, GL_STATIC_DRAW );
	int trinum = 0;
	for (int fi = 0; fi < p->faces; ++fi)
	{
	    face *f = &p->f[fi];
	    for (int i = 0; i < f->verts - 2; ++i)
	    {
		GLuint idx[3];
		idx[0] = f->vi[0];
		idx[1] = f->vi[i + 1];
		idx[2] = f->vi[i + 2];
		glBufferSubData( GL_ELEMENT_ARRAY_BUFFER, trinum++ * sizeof(idx), sizeof(idx), idx );
	    }
	}
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0 );
    }
    void bind( glsl::myprog *p )
    {
	glBindBuffer( GL_ARRAY_BUFFER, vbo[VTX] );
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, vbo[IDX] );

	for (int i = 0; i < 2; ++i)
	{
	    p->objvtx[i].bind( GL_FLOAT, sizeof(vtx), VBO_VTX_OFFSET_ARRAY(p,i) );
	    p->objvtx[i].enable();
	}
	p->objtex.bind( GL_FLOAT, sizeof(vtx), VBO_VTX_OFFSET(t) ); p->objtex.enable();
    }
    void unbind( glsl::myprog *p )
    {
	for (int i = 0; i < 2; ++i)
	    p->objvtx[i].disable();
	p->objtex.disable();

	glBindBuffer( GL_ARRAY_BUFFER, 0 );
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0 );
    }
    void draw()
    {
	glDrawElements( GL_TRIANGLES, 3 * tris, GL_UNSIGNED_INT, 0 );
    }
};

//
// globals
//

int nWinWidth = 1024;
int nWinHeight = 1024;
float fWinAspect = 1.0f;

bool bVerbose = false;
bool bAnimate = true;
bool bOutline = false;
bool bSwitchPolys = true;
bool bRandom = false;
bool bStereo = false;
#define MK_POLYHEDRON_ENUM(x) x,
enum Polyhedron { LIST_POLYHEDRA(MK_POLYHEDRON_ENUM) NUM_POLYHEDRONS } polyMode = Tetrahedron;
enum CHANGEMODE { ADD_PYRAMIDS, ADD_STELLATIONS } changeMode = ADD_PYRAMIDS;
enum KEEPMODE { KEEP_ALL, KEEP_POINTS, KEEP_SADDLES, KEEP_STELLATIONS } keepMode = KEEP_ALL;
int simplifyMode = -1;

float stellate_hz = 5;		// stellation frequency
float stellate = 0;		// stellation range (0..1)
float ax, ay, az;       	// angles for animation

poly		*polyTable[NUM_POLYHEDRONS][NUM_POLYHEDRONS] = { NULL };
GLuint		lastPolyVisit[NUM_POLYHEDRONS][NUM_POLYHEDRONS] = { 0 };
int		polyTableIndex = 0;
poly		*curpoly = NULL;
drawobj		*curobj = NULL;

// from shaders.glsl
glsl::myprog	*prog = NULL;

void selectPoly( Polyhedron frompoly, Polyhedron topoly );

//
// polyName - return a string name for the polyhedron type
//

const char *polyName( Polyhedron p )
{
    switch (p)
    {
    #define MK_POLYHEDRON_NAME(x) case x: return( #x );
    LIST_POLYHEDRA(MK_POLYHEDRON_NAME)
    }
    return( "(unknown)" );
}

//
// elapsed_ms - return milliseconds of elapsed time since last call
//

float elapsed_ms()
{
    static struct timeval otp;
    static bool bCalled = false;
    struct timeval tp;

    gettimeofday( &tp, NULL );

    float ms = bCalled ? (((tp.tv_sec - otp.tv_sec) * 1.0e3f) + ((tp.tv_usec - otp.tv_usec) * 1.0e-3f)) : 0.0f;
    bCalled = true;
    otp = tp;
    return( ms );
}

//
// randomizeByTime - use the current usec to seed the random generator
//

void randomizeByTime()
{
    struct timeval tp;
    gettimeofday( &tp, NULL );
    srandom( tp.tv_usec );
}

//
// animation - move the object around
//

void animation( float ms )
{
    ax += 0.5 / 16 * ms;
    ay += 0.3 / 16 * ms;
    az += 0.4 / 16 * ms;
    if (ax >= 360) ax = 0.0;
    if (ay >= 360) ay = 0.0;
    if (az >= 360) az = 0.0;

    stellate += stellate_hz * ms * 1e-3;
    if (bSwitchPolys)
    {
	if (stellate >= M_PI)
	{
	    stellate = fmodf( stellate, M_PI );
	    Polyhedron frompoly = Polyhedron(polyTableIndex / NUM_POLYHEDRONS);
	    Polyhedron topoly = Polyhedron(polyTableIndex % NUM_POLYHEDRONS);
	    ++lastPolyVisit[frompoly][topoly];
	    int polyVisit = ++lastPolyVisit[topoly][frompoly];
	    frompoly = topoly;
	    if (bRandom)
	    {
		int valid = 0;
		for (int i = 0; i < NUM_POLYHEDRONS; ++i)
		    valid += (lastPolyVisit[frompoly][i] < polyVisit) && (polyTable[frompoly][i] != NULL);
		int pick = random() % valid;
		for (int i = 0; i < NUM_POLYHEDRONS; ++i)
		{
		    if (polyTable[frompoly][i] && (lastPolyVisit[frompoly][i] < polyVisit) && (pick-- <= 0))
		    {
			topoly = Polyhedron(i);
			break;
		    }
		}
	    }
	    else // simple LRU
	    {
		GLuint oldest = polyVisit;
		for (int i = 0; i < NUM_POLYHEDRONS; ++i)
		{
		    if (polyTable[frompoly][i] && (lastPolyVisit[frompoly][i] < oldest))
		    {
			topoly = Polyhedron(i);
			oldest = lastPolyVisit[frompoly][i];
		    }
		}
	    }
	    selectPoly( frompoly, topoly );
	}
    }
    else
    {
	if (stellate >= 2 * M_PI)
	    stellate -= 2 * M_PI;
    }
}

//
// printMatrix - debugging aid
//

void printMatrix( const char *name, GLfloat m[] )
{
    printf( "%s (%f,%f,%f,%f),(%f,%f,%f,%f),(%f,%f,%f,%f),(%f,%f,%f,%f)\n", name,
	m[ 0], m[ 1], m[ 2], m[ 3],
	m[ 4], m[ 5], m[ 6], m[ 7],
	m[ 8], m[ 9], m[10], m[11],
	m[12], m[13], m[14], m[15] );
}

//
// drawString - render a string to the screen
//

void drawString( void *font, int x, int y, const char *text )
{
    unsigned char buff[1024];
    const char *s = text;
    unsigned char *d = buff;
    *d++ = *s++;
    while (*s)
    {
	if (isupper( *s )) *d++ = ' ';
	*d++ = *s++;
    }
    *d = '\0';

    y -= glutBitmapHeight( font ) / 2;
    x -= glutBitmapLength( font, buff ) / 2;
    glRasterPos2i( 0, 0 );
    glBitmap( 0, 0, 0, 0, x, y, NULL );
    glutBitmapString( font, buff );
}

//
// drawview - draw the object
//

void drawview( GLfloat blendfactor, GLfloat fEyeSeparation = 0.0f )
{
    prog->use();

    glClear( GL_DEPTH_BUFFER_BIT );

    glEnable( GL_DEPTH_TEST );

    GLfloat m[16];

    // HACK HACK HACK should use own matrix library...
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluPerspective( 10.0f, fWinAspect, 1.0f, 100.0f );

    glGetFloatv( GL_PROJECTION_MATRIX, m );
    prog->p.set( m );
    //printMatrix( "p ", m );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    gluLookAt( fEyeSeparation, 0, -20, 0, 0, 0, 0, 1, 0 );
    prog->eyepos.set( 0, 0, -20 );

    glRotatef( ax, 1.0, 0.0, 0.0 );
    glRotatef( ay, 0.0, 1.0, 0.0 );
    glRotatef( az, 0.0, 0.0, 1.0 );

    glGetFloatv( GL_MODELVIEW_MATRIX, m );
    prog->mv.set( m );
    //printMatrix( "mv ", m );

    curobj->bind( prog );
    curobj->draw();
    curobj->unbind( prog );

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glOrtho( 0, nWinWidth, 0, nWinHeight, -1, 1 );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glUseProgram( 0 );

    glDisable( GL_DEPTH_TEST );
    Polyhedron frompoly = Polyhedron(polyTableIndex / NUM_POLYHEDRONS);
    Polyhedron topoly = Polyhedron(polyTableIndex % NUM_POLYHEDRONS);
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glColor4f( 1, 1, 1, 1.0f - blendfactor );
    drawString( GLUT_BITMAP_HELVETICA_18, nWinWidth / 2, nWinHeight / 20, polyName( frompoly ) );
    glColor4f( 1, 1, 1, blendfactor );
    drawString( GLUT_BITMAP_HELVETICA_18, nWinWidth / 2, nWinHeight / 20, polyName( topoly ) );
    glDisable( GL_BLEND );
}

//
// drawscene - paint the screen
//

void drawscene()
{
    float ms = elapsed_ms();
    if (bAnimate) animation( ms );

    glViewport( 0, 0, nWinWidth, nWinHeight );

    if (bStereo) glClearColor( 0.1, 0.1, 0.1, 1.0 );
    else glClearColor( 0.2, 0.3, 0.1, 1.0 );
    glClear( GL_COLOR_BUFFER_BIT );

    glEnable( GL_CULL_FACE );
    //glDisable( GL_CULL_FACE );
    glCullFace( GL_BACK );

    glDepthFunc( GL_LEQUAL );

    prog->use();
    prog->outline.set( bOutline );
    prog->luma.set( bStereo );

    GLfloat blendfactor = 0.5f * (1.0f - cosf( stellate ));
    prog->stellate.set( blendfactor );

    if (bStereo)
    {
	glColorMask( GL_TRUE, GL_FALSE, GL_FALSE, GL_FALSE );
	drawview( blendfactor, 1.0f );
	glColorMask( GL_FALSE, GL_TRUE, GL_TRUE, GL_FALSE );
	drawview( blendfactor, -1.0f );
	glColorMask( GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE );
    }
    else drawview( blendfactor );

    glutSwapBuffers();
    glutPostRedisplay();
}

//
// changePoly - generated a stellated version of the shape
//

void changePoly( int tovtx, float scale )
{
    poly *newpoly = NULL;
    switch (changeMode)
    {
    case ADD_PYRAMIDS:
	newpoly = addPyramids( curpoly, tovtx, scale );
	break;
    case ADD_STELLATIONS:
	newpoly = addStellations( curpoly, tovtx, scale );
	break;
    }
    if (newpoly && (keepMode != KEEP_ALL))
    {
	poly *keeper = NULL;
	switch (keepMode)
	{
	case KEEP_POINTS:
	    keeper = keepVertexType( newpoly, tovtx ^ 1, vtx::POINT );
	    break;
	case KEEP_SADDLES:
	    keeper = keepVertexType( newpoly, tovtx ^ 1, vtx::SADDLE );
	    break;
	case KEEP_STELLATIONS:
	    keeper = keepVertexType( newpoly, tovtx ^ 1, vtx::STELLATE );
	    break;
	}
	if (keeper)
	{
	    delete newpoly;
	    newpoly = keeper;
	    keeper = NULL;
	}
    }
    if (newpoly)
    {
	delete curpoly;
	curpoly = newpoly;
    }
    if (simplifyMode >= 0)
    {
	#if defined(VERBOSE_SIMPLIFY)
	    printf( "simplifying %d verts, %d faces:\n    ", curpoly->verts, curpoly->next_face );
	    for (int i = 0; i < curpoly->verts; ++i)
		printf( " %2d", i );
	    printf( "\n" );
	    for (int i = 0; i < curpoly->verts; ++i)
	    {
		printf( " %2d:", i );
		for (int j = 0; j < curpoly->verts; ++j)
		{
		    int fi = curpoly->faceidx( i, j );
		    if (fi >= 0) printf( " %2d", fi ); else printf( "  ." );
		}
		printf( "\n" );
	    }
	#endif
	newpoly = simplify( curpoly, simplifyMode );
	delete curpoly;
	curpoly = newpoly;
    }
    curobj->set( curpoly );
}

//
// makeTetrahedron - generate a poly object for a tetrahedron
//

poly *makeTetrahedron()
{
    static const float a = powf( 1.0f / 3.0f, 0.5f );

    static vtx v[4] =
    {
	{ { {  a,  a,  a }, {  a,  a,  a } }, { 1, 1 }, vtx::POINT },
	{ { {  a, -a, -a }, {  a, -a, -a } }, { 1, 0 }, vtx::POINT },
	{ { { -a, -a,  a }, { -a, -a,  a } }, { 0, 1 }, vtx::POINT },
	{ { { -a,  a, -a }, { -a,  a, -a } }, { 0, 0 }, vtx::POINT },
    };

    static int idx[4][3] =
    {
	{ 0, 2, 1 },
	{ 2, 3, 1 },
	{ 3, 0, 1 },
	{ 0, 3, 2 }
    };

    poly *p = new poly( 4, 4 );
    for (int i = 0; i < 4; ++i) p->setvtx( i, v[i] );
    for (int i = 0; i < 4; ++i) p->addface( 3, idx[i] );
    return( p );
}

//
// setupTetrahedron - generate poly and drawobj objects for a tetrahedron
//

void setupTetrahedron()
{
    delete curpoly;
    curpoly = makeTetrahedron();
    changePoly( 1, 3.0f );
}

//
// setupCube - generate poly and draw objects for a cube
//

void setupCube()
{
    static const float a = powf( 1.0f / 3.0f, 0.5f );

    static vtx v[8] =
    {
	{ { {  a,  a,  a }, {  a,  a,  a } }, { 1, 1 }, vtx::POINT },
	{ { {  a, -a,  a }, {  a, -a,  a } }, { 1, 0 }, vtx::POINT },
	{ { {  a, -a, -a }, {  a, -a, -a } }, { 0, 0 }, vtx::POINT },
	{ { {  a,  a, -a }, {  a,  a, -a } }, { 0, 1 }, vtx::POINT },
	{ { { -a,  a,  a }, { -a,  a,  a } }, { 1, 1 }, vtx::POINT },
	{ { { -a, -a,  a }, { -a, -a,  a } }, { 1, 0 }, vtx::POINT },
	{ { { -a, -a, -a }, { -a, -a, -a } }, { 0, 0 }, vtx::POINT },
	{ { { -a,  a, -a }, { -a,  a, -a } }, { 0, 1 }, vtx::POINT },
    };

    static int idx[6][4] =
    {
	{ 0, 1, 2, 3 },
	{ 0, 4, 5, 1 },
	{ 4, 7, 6, 5 },
	{ 2, 6, 7, 3 },
	{ 0, 3, 7, 4 },
	{ 1, 5, 6, 2 }
    };

    delete curpoly;
    curpoly = new poly( 8, 6 );
    for (int i = 0; i < 8; ++i) curpoly->setvtx( i, v[i] );
    for (int i = 0; i < 6; ++i) curpoly->addface( 4, idx[i] );

    changePoly( 0, 2.0f );
}

//
// setPolyMode - select the base polygon
//

void setPolyMode( Polyhedron p )
{
    polyMode = p;
    switch (p)
    {
    case Tetrahedron:
	setupTetrahedron();
	break;
    case Cube:
	setupCube();
	break;
    }
}

//
// setChangeMode - set the morph type
//

void setChangeMode( CHANGEMODE c )
{
    changeMode = c;
    setPolyMode( polyMode );
}

//
// setKeepMode - set the morph type
//

void setKeepMode( KEEPMODE k )
{
    keepMode = k;
    setPolyMode( polyMode );
}

//
// toggleSimplify - switch simplification on/off
//

void toggleSimplify()
{
    if (++simplifyMode > 1) simplifyMode = -1;
    setPolyMode( polyMode );
}

//
// selectPoly - advance through the polyTable for debugging
//

void selectPoly( int dir )
{
    poly *p = NULL;
    for (int i = 0; i < NUM_POLYHEDRONS * NUM_POLYHEDRONS; ++i)
    {
	polyTableIndex += dir;
	dir += (dir == 0);
	if (polyTableIndex < 0) polyTableIndex = NUM_POLYHEDRONS * NUM_POLYHEDRONS - 1;
	else if (polyTableIndex >= NUM_POLYHEDRONS * NUM_POLYHEDRONS) polyTableIndex = 0;
	Polyhedron frompoly = Polyhedron(polyTableIndex / NUM_POLYHEDRONS);
	Polyhedron topoly = Polyhedron(polyTableIndex % NUM_POLYHEDRONS);
	p = polyTable[frompoly][topoly];
	if (p)
	{
	    curobj->set( p );
	    if (bVerbose)
	    {
		if (frompoly == topoly) printf( "%s: ", polyName( frompoly ) );
		else printf( "%s->%s: ", polyName( frompoly ), polyName( topoly ) );
		printf( " %d verts, %d faces, %d edges\n", p->verts, p->next_face, p->edges() );
	    }
	    break;
	}
    }
}

void selectPoly( Polyhedron frompoly, Polyhedron topoly )
{
    polyTableIndex = frompoly * NUM_POLYHEDRONS + topoly;
    selectPoly( 0 );
}

//
// generatePolyTable - fill a table of transitions from one polyhedron to the others
//

void generatePolyTable()
{
    #if defined(VERBOSE_SIMPLIFY)
	#define MK_POLY(f,t) \
	printf( "Making %s->%s\n", polyName( f ), polyName( t ) ); \
	polyTable[f][t]
    #else
	#define MK_POLY(f,t) polyTable[f][t]
    #endif

    MK_POLY(Tetrahedron,Tetrahedron) = makeTetrahedron();
    MK_POLY(Tetrahedron,StellaOctangula) = addStellations( polyTable[Tetrahedron][Tetrahedron], 1 );
    MK_POLY(Tetrahedron,Cube) = addPyramids( polyTable[Tetrahedron][Tetrahedron], 1 );

    MK_POLY(StellaOctangula,StellaOctangula) = simplify( polyTable[Tetrahedron][StellaOctangula], 1 );
    MK_POLY(StellaOctangula,Octahedron) = keepVertexType( polyTable[StellaOctangula][StellaOctangula], 1, vtx::SADDLE );

    MK_POLY(Cube,Cube) = simplify( polyTable[Tetrahedron][Cube], 1, true );
    MK_POLY(Cube,Dodecahedron) = addTents( polyTable[Cube][Cube], 1, 0.0f, M_PHI );
    MK_POLY(Cube,CubeOctahedronCompound) = addStellations( polyTable[Cube][Cube], 1 );
    MK_POLY(Cube,RhombicDodecahedron) = addPyramids( polyTable[Cube][Cube], 1 );

    MK_POLY(CubeOctahedronCompound,CubeOctahedronCompound) = simplify( polyTable[Cube][CubeOctahedronCompound], 1 );
    MK_POLY(CubeOctahedronCompound,Cuboctahedron) = keepVertexType( polyTable[CubeOctahedronCompound][CubeOctahedronCompound], 1, vtx::SADDLE );

    MK_POLY(Octahedron,Octahedron) = simplify( polyTable[StellaOctangula][Octahedron], 1, true );
    MK_POLY(Octahedron,CubeOctahedronCompound) = addStellations( polyTable[Octahedron][Octahedron], 1 );
    MK_POLY(Octahedron,RhombicDodecahedron) = addPyramids( polyTable[Octahedron][Octahedron], 1 );
    MK_POLY(Octahedron,Icosahedron) = removeTents( polyTable[Octahedron][Octahedron], 1 );

    MK_POLY(RhombicDodecahedron,RhombicDodecahedron) = simplify( polyTable[Cube][RhombicDodecahedron], 1, true );
    MK_POLY(RhombicDodecahedron,DeltoidalIcositetrahedron) = addPyramids( polyTable[RhombicDodecahedron][RhombicDodecahedron], 1 );

    MK_POLY(Cuboctahedron,Cuboctahedron) = simplify( polyTable[CubeOctahedronCompound][Cuboctahedron], 1, true );
    MK_POLY(Cuboctahedron,DeltoidalIcositetrahedron) = addPyramids( polyTable[Cuboctahedron][Cuboctahedron], 1 );

    MK_POLY(Dodecahedron,Dodecahedron) = simplify( polyTable[Cube][Dodecahedron], 1, true );
    MK_POLY(Dodecahedron,DodecahedronIcosahedronCompound) = addStellations( polyTable[Dodecahedron][Dodecahedron], 1 );
    MK_POLY(Dodecahedron,RhombicDodecahedron) = addTents( polyTable[Cube][Cube], 1, M_PHI, 1.0f );
    MK_POLY(Dodecahedron,RhombicTriacontahedron) = addPyramids( polyTable[Dodecahedron][Dodecahedron], 1 );

    MK_POLY(DodecahedronIcosahedronCompound,DodecahedronIcosahedronCompound) = simplify( polyTable[Dodecahedron][DodecahedronIcosahedronCompound], 1 );
    MK_POLY(DodecahedronIcosahedronCompound,Icosahedron) = keepVertexType( polyTable[DodecahedronIcosahedronCompound][DodecahedronIcosahedronCompound], 1, vtx::STELLATE );
    MK_POLY(DodecahedronIcosahedronCompound,Icosidodecahedron) = keepVertexType( polyTable[DodecahedronIcosahedronCompound][DodecahedronIcosahedronCompound], 1, vtx::SADDLE );

    MK_POLY(Icosahedron,Icosahedron) = simplify( polyTable[DodecahedronIcosahedronCompound][Icosahedron], 1, true );
    MK_POLY(Icosahedron,RhombicTriacontahedron) = addPyramids( polyTable[Icosahedron][Icosahedron], 1 );

    MK_POLY(RhombicTriacontahedron,RhombicTriacontahedron) = simplify( polyTable[Dodecahedron][RhombicTriacontahedron], 1, true );
    MK_POLY(RhombicTriacontahedron,DeltoidalHexecontahedron) = addPyramids( polyTable[RhombicTriacontahedron][RhombicTriacontahedron], 1 );

    MK_POLY(Icosidodecahedron,Icosidodecahedron) = simplify( polyTable[DodecahedronIcosahedronCompound][Icosidodecahedron], 1, true );
    MK_POLY(Icosidodecahedron,DeltoidalHexecontahedron) = addPyramids( polyTable[Icosidodecahedron][Icosidodecahedron], 1 );

    MK_POLY(DeltoidalIcositetrahedron,DeltoidalIcositetrahedron) = simplify( polyTable[RhombicDodecahedron][DeltoidalIcositetrahedron], 1 );

    MK_POLY(DeltoidalHexecontahedron,DeltoidalHexecontahedron) = simplify( polyTable[RhombicTriacontahedron][DeltoidalHexecontahedron], 1 );

    for (int i = 0; i < NUM_POLYHEDRONS; ++i)
    {
	for (int j = 0; j < NUM_POLYHEDRONS; ++j)
	{
	    if (polyTable[i][j] && !polyTable[j][i])
		polyTable[j][i] = flip( polyTable[i][j] );
	}
    }

    selectPoly( Tetrahedron, Tetrahedron );

    if (bVerbose)
    {
	printf( "polyTable:\n" );
	for (int i = 0; i < NUM_POLYHEDRONS; ++i)
	{
	    printf( "%35s:", polyName( Polyhedron(i) ) );
	    for (int j = 0; j < NUM_POLYHEDRONS; ++j)
		printf( polyTable[i][j] ? " *" : " ." );
	    printf( "\n" );
	}
    }
}

//
// initLightAndMaterial - init everything
//

void initLightAndMaterial()
{
    prog = new glsl::myprog();
    char buff[10240];
    GLsizei len = sizeof(buff);
    if (!prog->build( buff, &len ))
    {
	WARN(( "Failed to build myprog:\n%s\n", buff ));
    }
    if (bVerbose)
    {
	len = sizeof(buff); prog->myvtx.info( buff, &len );
	printf( "myvtx messages:\n%s\n", buff );
	len = sizeof(buff); prog->mygeom.info( buff, &len );
	printf( "mygeom messages:\n%s\n", buff );
	len = sizeof(buff); prog->myfrag.info( buff, &len );
	printf( "myfrag messages:\n%s\n", buff );
	len = sizeof(buff); prog->info( buff, &len );
	printf( "linker messages:\n%s\n", buff );
    }

    prog->lightpos.set( 10, 10, -5 );
    prog->lightcolor.set( 1, 1, 1 );
    prog->objcolor.set( 1, 0, 0 );
    prog->warmcolor.set( 1, 1, 0 );
    prog->coolcolor.set( 0, 0, 1 );
    prog->ambient.set( 0.1f );
}

//
// resize - handle window resize events
//

void resize( int w, int h )
{
    nWinWidth = w;
    nWinHeight = h;
    fWinAspect = float(w) / (h ? h : 1);
}

//
// command - handle keyboard, menu or special key commands
//

void command( int cmd )
{
    switch (cmd)
    {
    #define MK_CMD_CASE(label,value,case,cmd) case cmd; break;
    LIST_COMMANDS(MK_CMD_CASE)
    }
}

//
// keyboard - handle GLUT keypresses
//

void keyboard( unsigned char c, int x, int y )
{
    command( c );
}

//
// special - handle other GLUT keypresses
//

void special( int c, int x, int y )
{
    command( MK_SPECIALKEY(c) );
}

//
// main - make it all happen
//

int main( int argc, char *argv[] )
{
    glutInit( &argc, argv );

    glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE );
    glutInitWindowSize( nWinWidth, nWinHeight );
    glutCreateWindow( "Platonic" );

    ax = 0.0;
    ay = 0.0;
    az = 0.0;

    randomizeByTime();

    initLightAndMaterial();

    curobj = new drawobj();

    generatePolyTable();

    glutDisplayFunc( drawscene );
    glutReshapeFunc( resize );
    glutCreateMenu( command );
    #define MK_MENU(label,value,case,cmd) glutAddMenuEntry( label, value );
    LIST_COMMANDS(MK_MENU)
    glutAttachMenu( GLUT_RIGHT_BUTTON );
    glutKeyboardFunc( keyboard );
    glutSpecialFunc( special );

    glutMainLoop();

    return( 0 );
}
