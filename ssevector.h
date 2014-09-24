//
// ssevector.h - implement float vector types via SSE intrinsics
//

#if !defined(SSEVECTOR_H)
#define SSEVECTOR_H

#include <math.h>
#include <xmmintrin.h>
#include <emmintrin.h>

#define ALIGNSSE __declspec(align(16))
static const int SSE_FLOAT_SIZE = 4;
static const int SSE_FLOAT_MASK = SSE_FLOAT_SIZE - 1;

#if !defined(M_PI)
    #define M_PI 3.14159265358979323846
#endif

//
// svec4 - memory-aligned SSE 4-vector
//

class svec4
{
protected:
    __m128	r;
    operator __m128&()
    {
	return this->r;
    }
    svec4( __m128 x ) : r(x) {}
    operator float*() const
    {
	return static_cast<float *>(static_cast<void *>(const_cast<__m128 *>(&this->r)));
    }

public:
    svec4() { r = _mm_setzero_ps(); }
    svec4( float x ) { r = _mm_set_ps1( x ); }
    svec4( float x, float y, float z, float w ) { r = _mm_set_ps( w, z, y, x ); }
    svec4( float *p ) { r = _mm_load_ps( p ); }
    svec4( const svec4 &v ) : r(v.r) {}
    svec4 &operator=(const svec4 &v)
    {
	_mm_store_ps( (float *)(this), svec4(v) );
	return *this;
    }
    void store(float *p) const
    {
	_mm_store_ps( p, svec4(*this) );
    }
    void storeu(float *p) const
    {
	_mm_storeu_ps( p, svec4(*this) );
    }
    float &operator[]( int i ) const
    {
	return ((float *)(this))[i];
    }
    operator float() const
    {
	return ((float *)(this))[0];
    }
    // arithmetic operators
    svec4 abs() const
    {
	return svec4(_mm_max_ps( svec4(*this), -svec4(*this) ));
	//return select( *this >= svec4(), -svec4(*this) );
    }
    svec4 rcp_approx() const
    {
	return svec4(_mm_rcp_ps( svec4(*this) ));
    }
    svec4 rcp() const
    {
	// use newton's method to refine rcp
	svec4 x = rcp_approx();
	return x * (svec4(2.0f) - svec4(*this) * x);
    }
    float rcp_approx_ss() const
    {
	return svec4(_mm_rcp_ss( svec4(*this) ))[0];
    }
    float rcp_ss() const
    {
	// use newton's method to refine rcp
	svec4 x( rcp_approx_ss() );
	return svec4(_mm_mul_ss( x, _mm_sub_ss( svec4(2.0f), _mm_mul_ss( svec4(*this), x ) ) ))[0];
    }
    svec4 sqrt() const
    {
	return svec4(_mm_sqrt_ps( svec4(*this) ));
    }
    float sqrt_ss() const
    {
	return svec4(_mm_sqrt_ss( svec4(*this) ))[0];
    }
    svec4 rsqrt_approx() const
    {
	return svec4(_mm_rsqrt_ps( svec4(*this) ));
    }
    svec4 rsqrt() const
    {
	// use newton's method to refine rsqrt
	svec4 x = rsqrt_approx();
	return svec4(-0.5f) * x * (x * x * svec4(*this) - svec4(3.0f));
    }
    float rsqrt_approx_ss() const
    {
	return svec4(_mm_rsqrt_ss( svec4(*this) ))[0];
    }
    float rsqrt_ss() const
    {
	// use newton's method to refine rsqrt
	svec4 x( rsqrt_approx_ss() );
	svec4 x2r_3( _mm_sub_ss( _mm_mul_ss( _mm_mul_ss( x, x ), svec4(*this) ), svec4(3.0f) ) );
	return svec4(_mm_mul_ss( _mm_mul_ss( svec4(-0.5f), x ), x2r_3 ))[0];
    }
    float length() const
    {
	return dot( *this ).sqrt_ss();
    }
    svec4 operator -() const
    {
	return svec4(_mm_sub_ps( svec4(), svec4(*this) ));
    }
    svec4 dot(const svec4 &v) const
    {
	svec4 t = *this * v;
	// 3 2 1 0 + 3 3 1 1 = 3+3 3+2 1+1 1+0
	svec4 t2 = _mm_add_ps( _mm_shuffle_ps( t, t, _MM_SHUFFLE(3,3,1,1) ), t );
	// 3+3 3+2 1+1 1+0 + 3+2 3+2 3+2 3+2 = x x x 3+2+1+0
	svec4 t3 = _mm_add_ss( _mm_shuffle_ps( t2, t2, _MM_SHUFFLE(2,2,2,2) ), t2 );
	return svec4(t3[0]);
    }
    svec4 operator + (const svec4 &v) const
    {
	return svec4(_mm_add_ps( svec4(*this), svec4(v) ));
    }
    svec4 operator - (const svec4 &v) const
    {
	return svec4(_mm_sub_ps( svec4(*this), svec4(v) ));
    }
    svec4 operator * (const svec4 &v) const
    {
	return svec4(_mm_mul_ps( svec4(*this), svec4(v) ));
    }
    svec4 operator / (const svec4 &v) const
    {
	return svec4(_mm_div_ps( svec4(*this), svec4(v) ));
    }
    // operators which modify this value
    svec4 &operator += (const svec4 &v)
    {
	return *this = *this + v;
    }
    svec4 &operator -= (const svec4 &v)
    {
	return *this = *this - v;
    }
    svec4 &operator *= (const svec4 &v)
    {
	return *this = *this * v;
    }
    svec4 &operator /= (const svec4 &v)
    {
	return *this = *this / v;
    }
    svec4 &normalize()
    {
	float n = dot( *this ).rsqrt_ss();
	return *this *= n;
    }
    svec4 &clamp()
    {
	return *this = max( svec4() );
    }
    svec4 &saturate()
    {
	return *this = max( svec4() ).min( svec4( 1.0f ) );
    }
    // logical operators
    operator bool() const
    {
	// return false if all elements are zero, true otherwise
	return bool(_mm_movemask_ps( r ) != 0);
    }
    svec4 operator ! () const
    {
	return svec4(_mm_cmpeq_ps( svec4(*this), svec4() ));
    }
    svec4 operator == (const svec4 &v) const
    {
	return svec4(_mm_cmpeq_ps( svec4(*this), svec4(v) ));
    }
    svec4 operator != (const svec4 &v) const
    {
	return svec4(_mm_cmpneq_ps( svec4(*this), svec4(v) ));
    }
    svec4 operator < (const svec4 &v) const
    {
	return svec4(_mm_cmplt_ps( svec4(*this), svec4(v) ));
    }
    svec4 operator <= (const svec4 &v) const
    {
	return svec4(_mm_cmple_ps( svec4(*this), svec4(v) ));
    }
    svec4 operator > (const svec4 &v) const
    {
	return svec4(_mm_cmpgt_ps( svec4(*this), svec4(v) ));
    }
    svec4 operator >= (const svec4 &v) const
    {
	return svec4(_mm_cmpge_ps( svec4(*this), svec4(v) ));
    }
    svec4 min(const svec4 &v) const
    {
	return svec4(_mm_min_ps( svec4(*this), svec4(v) ));
    }
    svec4 max(const svec4 &v) const
    {
	return svec4(_mm_max_ps( svec4(*this), svec4(v) ));
    }
    // masking
    svec4 operator ~ () const
    {
	// HACK HACK HACK what's the fastest/safest way to get all 1's? - Ed
	return svec4(_mm_andnot_ps( svec4(*this), _mm_cmpeq_ps( svec4(), svec4() ) ));
    }
    svec4 operator & (const svec4 &v) const
    {
	return svec4(_mm_and_ps( svec4(*this), svec4(v) ));
    }
    svec4 operator | (const svec4 &v) const
    {
	return svec4(_mm_or_ps( svec4(*this), svec4(v) ));
    }
    svec4 operator ^ (const svec4 &v) const
    {
	return svec4(_mm_xor_ps( svec4(*this), svec4(v) ));
    }
    svec4 &operator &= (const svec4 &v)
    {
	return *this = *this & v;
    }
    svec4 &operator |= (const svec4 &v)
    {
	return *this = *this | v;
    }
    svec4 &operator ^= (const svec4 &v)
    {
	return *this = *this ^ v;
    }
    // select this where mask is true, v where mask is false
    svec4 select(const svec4 &mask, const svec4 &v) const
    {
	// note that andnot(a,b) is really notand: (~a & b)
	return svec4(_mm_or_ps( _mm_and_ps( svec4(mask), svec4(*this) ), _mm_andnot_ps( svec4(mask), svec4(v) ) ));
    }
    // accumulate into a memory location
    const svec4 &accumulate(float *p) const
    {
	// avoid memory transaction if the value is 0
	if (bool(*this != svec4()))
	{
	    svec4 t( p );
	    t += *this;
	    t.store( p );
	}
	return *this;
    }
};

#endif // defined(SSEVECTOR_H)
