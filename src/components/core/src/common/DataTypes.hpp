/*
 * DataTypes.hpp
 *
 *  Created on: Jun 16, 2016
 *      Author: rrsettgast
 */

#ifndef DATATYPES_HPP
#define DATATYPES_HPP

#include <cassert>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>
#include <typeindex>
#include <typeinfo>
#include <map>
#include <unordered_map>
#include <vector>


#include "Macros.hpp"


#include "legacy/ArrayT/ArrayT.h"
#include "math/TensorT/TensorT.h"

#ifndef CONTAINERARRAY_RETURN_PTR
#define CONTAINERARRAY_RETURN_PTR 1
#endif


namespace geosx
{

using  int32      = std::int32_t;
using uint32      = std::uint32_t;
using  int64      = std::int64_t;
using uint64      = std::uint64_t;
using std_size_t  = std::size_t;
using string      = std::string;
using integer     = int;
using uinteger    = unsigned int;
using index_t     = int;


using real32 = float;
using real64 = double;
using real   = double;

template< typename T >
using ptr = T*;

template< typename T >
using c_ptr = T const *;


using int32_ptr        = ptr<int32>;
using int32_const_ptr  = c_ptr<int32>;

using uint32_ptr        = ptr<uint32>;
using uint32_const_ptr  = c_ptr<uint32>;

using int64_ptr        = ptr<int64>;
using int64_const_ptr  = c_ptr<int64>;

using uint64_ptr        = ptr<uint64>;
using uint64_const_ptr  = c_ptr<uint64>;

using real32_ptr        = ptr<real32>;
using real32_const_ptr  = c_ptr<real32>;

using real64_ptr        = ptr<real64>;
using real64_const_ptr  = c_ptr<real64>;

//***** BEGIN LEGACY TYPEDEFS *****
using globalIndex = int64;
using localIndex  = int32;

using realT    = double;


template< typename T >
//using array = std::vector<T>;
using array = Array1dT<T>;

template< typename T >
using set = std::set<T>;

template< typename TKEY, typename TVAL >
using map = std::map<TKEY,TVAL>;

template< typename TKEY, typename TVAL >
using unordered_map = std::unordered_map<TKEY,TVAL>;


//***** END LEGACY TYPEDEFS *****



using int32_array        = array<int32>;
using int32_const_array  = array<int32 const>;

using uint32_array        = array<uint32>;
using uint32_const_array  = array<uint32 const>;

using int64_array        = array<int64>;
using int64_const_array  = array<int64 const>;

using uint64_array        = array<uint64>;
using uint64_const_array  = array<uint64 const>;

using real32_array        = array<real32>;
using real32_const_array  = array<real32 const>;

using real64_array        = array<real64>;
using real64_const_array  = array<real64 const>;

using string_array        = array<string>;
using string_const_array  = array<string const>;

using int32_set        = set<int32>;
using int32_const_set  = set<int32 const>;

using uint32_set        = set<uint32>;
using uint32_const_set  = set<uint32 const>;

using int64_set        = set<int64>;
using int64_const_set  = set<int64 const>;

using uint64_set        = set<uint64>;
using uint64_const_set  = set<uint64 const>;

using real32_set        = set<real32>;
using real32_const_set  = set<real32 const>;

using real64_set        = set<real64>;
using real64_const_set  = set<real64 const>;

using string_set        = set<string>;
using string_const_set  = set<string const>;



//***** BEGIN LEGACY TYPEDEFS *****
using rArray1d = Array1dT<real64>;
using iArray1d = Array1dT<int32>;
using lArray1d = Array1dT<localIndex>;
using gArray1d = Array1dT<globalIndex>;

typedef Array1dT<std::string> sArray1d;

typedef std::set<localIndex> lSet;
typedef std::set<globalIndex> gSet;

typedef int FieldKey;

typedef Array2dT<localIndex> lArray2d;
typedef Array1dT<std::pair<int,localIndex> > pArray1d;
typedef std::set<std::pair<int,localIndex> > pSet;

using r1_array = array<R1Tensor>;
using r2_array = array<R2Tensor>;
using r2Sym_array = array<R2SymTensor>;

//using mapPair = std::pair<int32, localIndex>;
using mapPair_array = std::pair<int32_array, int32_array>;
//***** END LEGACY TYPEDEFS *****

class rtTypes
{
public:

  static std::string typeNames( std::type_index const key )
  {
    const std::unordered_map<std::type_index, std::string> type_names =
    {
      {std::type_index(typeid(int32)), "int32"},
      {std::type_index(typeid(uint32)), "uint32"},
      {std::type_index(typeid(int64)), "int64"},
      {std::type_index(typeid(uint64)), "uint64"},
      {std::type_index(typeid(real32)), "real32"},
      {std::type_index(typeid(real64)), "real64"},
      {std::type_index(typeid(R1Tensor)), "r1Tensor"},
      {std::type_index(typeid(R2Tensor)), "r2Tensor"},
      {std::type_index(typeid(R2SymTensor)), "r2SymTensor"},
      {std::type_index(typeid(int32_array)), "int32_array"},
      {std::type_index(typeid(uint32_array)), "uint32_array"},
      {std::type_index(typeid(int64_array)), "int64_array"},
      {std::type_index(typeid(uint64_array)), "uint64_array"},
      {std::type_index(typeid(real32_array)), "real32_array"},
      {std::type_index(typeid(real64_array)), "real64_array"},
      {std::type_index(typeid(r1_array)), "r1_array"},
      {std::type_index(typeid(r2_array)), "r2_array"},
      {std::type_index(typeid(r2Sym_array)), "r2Sym_array"},
      {std::type_index(typeid(std_size_t)), "std_size_t"},
      {std::type_index(typeid(string)), "string"},
      {std::type_index(typeid(mapPair_array)), "mapPair_array"}
    };
    return type_names.at(key);
  }



  enum class TypeIDs
  {
    int32_id,
    uint32_id,
    int64_id,
    uint64_id,
    real32_id,
    real64_id,
    r1Tensor_id,
    r2Tensor_id,
    r2SymTensor_id,
    int32_array_id,
    uint32_array_id,
    int64_array_id,
    uint64_array_id,
    real32_array_id,
    real64_array_id,
    r1_array_id,
    r2_array_id,
    r2Sym_array_id,
    std_size_t_id,
    string_id,
    string_array_id,
    mapPair_array_id,
    none_id
  };

  static TypeIDs typeID( string const & name )
  {
    const std::unordered_map<string,TypeIDs> type_names =
    {
      { "int32",        TypeIDs::int32_id },
      { "uint32",       TypeIDs::uint32_id },
      { "int64",        TypeIDs::int64_id },
      { "uint64",       TypeIDs::uint64_id },
      { "real32",       TypeIDs::real32_id },
      { "real64",       TypeIDs::real64_id },
      { "R1Tensor",     TypeIDs::r1Tensor_id },
      { "R2Tensor",     TypeIDs::r2Tensor_id },
      { "R2SymTensor",  TypeIDs::r2SymTensor_id },
      { "int32_array",  TypeIDs::int32_array_id },
      { "uint32_array", TypeIDs::uint32_array_id },
      { "int64_array",  TypeIDs::int64_array_id },
      { "uint64_array", TypeIDs::uint64_array_id },
      { "real32_array", TypeIDs::real32_array_id },
      { "real64_array", TypeIDs::real64_array_id },
      { "r1_array",     TypeIDs::r1_array_id },
      { "r2_array",     TypeIDs::r2_array_id },
      { "r2Sym_array",  TypeIDs::r2Sym_array_id },
      { "std_size_t",   TypeIDs::std_size_t_id },
      { "string",       TypeIDs::string_id },
      { "string_array", TypeIDs::string_array_id },
      { "mapPair_array",      TypeIDs::mapPair_array_id },
      { "",             TypeIDs::none_id }
    };
    return type_names.at(name);
  }

  static TypeIDs typeID( std::type_index typeIndex )
  {
    const std::unordered_map<std::type_index,TypeIDs> type_names =
    {
      { std::type_index(typeid(int32)),        TypeIDs::int32_id },
      { std::type_index(typeid(uint32)),       TypeIDs::uint32_id },
      { std::type_index(typeid(int64)),        TypeIDs::int64_id },
      { std::type_index(typeid(uint64)),       TypeIDs::uint64_id },
      { std::type_index(typeid(real32)),       TypeIDs::real32_id },
      { std::type_index(typeid(real64)),       TypeIDs::real64_id },
      { std::type_index(typeid(R1Tensor)),     TypeIDs::r1Tensor_id },
      { std::type_index(typeid(R2Tensor)),     TypeIDs::r2Tensor_id },
      { std::type_index(typeid(R2SymTensor)),  TypeIDs::r2SymTensor_id },
      { std::type_index(typeid(int32_array)),  TypeIDs::int32_array_id },
      { std::type_index(typeid(uint32_array)), TypeIDs::uint32_array_id },
      { std::type_index(typeid(int64_array)),  TypeIDs::int64_array_id },
      { std::type_index(typeid(uint64_array)), TypeIDs::uint64_array_id },
      { std::type_index(typeid(real32_array)), TypeIDs::real32_array_id },
      { std::type_index(typeid(real64_array)), TypeIDs::real64_array_id },
      { std::type_index(typeid(r1_array)),     TypeIDs::r1_array_id },
      { std::type_index(typeid(r2_array)),     TypeIDs::r2_array_id },
      { std::type_index(typeid(r2Sym_array)),  TypeIDs::r2Sym_array_id },
      { std::type_index(typeid(std_size_t)),   TypeIDs::std_size_t_id },
      { std::type_index(typeid(string)),       TypeIDs::string_id },
      { std::type_index(typeid(string_array)), TypeIDs::string_array_id },
      { std::type_index(typeid(mapPair_array)),TypeIDs::mapPair_array_id }
    };
    return type_names.at(typeIndex);
  }



  // Matching regex for data types in xml
  class typeRegex
  {
  private:
    std::string ru = "[0-9]*";
    std::string ri = "[+-]?[0-9]*";    
    std::string rr = "[0-9]*\\.?([0-9]*)?[eE]?[-+]?([0-9]*)?";
    std::string rs = "[a-zA-Z0-9_,\\(\\)+-/\\*]*";
    std::string r1 = rr + ",? " + rr + ",? " + rr;
    std::string r2 = rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr;
    std::string r2s = rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr;
    
    std::unordered_map<std::string, std::string> regexMap = 
    {
      {"int32", ri},
      {"uint32", ru},
      {"int64", ri},
      {"uint64", ru},
      {"real32", rr},
      {"real64", rr},
      {"R1Tensor", r1},
      {"R2Tensor", r2},
      {"R2SymTensor", r2s},
      {"int32_array", "((" + ri + ",? )*)?" + ri},
      {"uint32_array", "((" + ru + ",? )*)?" + ru},
      {"int64_array", "((" + ri + ",? )*)?" + ri},
      {"uint64_array", "((" + ru + ",? )*)?" + ru},
      {"real32_array", "((" + rr + ",? )*)?" + rr},
      {"real64_array", "((" + rr + ",? )*)?" + rr},
      {"r1_array", "((" + r1 + "; )*)?" + r1},
      {"r2_array", "((" + r2 + "; )*)?" + r2},
      {"r2Sym_array", "((" + r2s + "; )*)?" + r2s},
      {"std_size_t", ru},
      {"string", rs},
      {"string_array", "((" + rs + ",? )*)?" + rs},
      {"mapPair", rs},
      {"mapPair_array", "((" + rs + ",? )*)?" + rs}
    };

  public:
    std::unordered_map<std::string, std::string>::iterator begin(){return regexMap.begin();}
    std::unordered_map<std::string, std::string>::iterator end(){return regexMap.end();}
    std::unordered_map<std::string, std::string>::const_iterator begin() const {return regexMap.begin();}
    std::unordered_map<std::string, std::string>::const_iterator end() const {return regexMap.end();}
  };



  template< typename LAMBDA >
  static auto ApplyIntrinsicTypeLambda1( const TypeIDs type,
                                         LAMBDA lambda )
  {
    switch( type )
    {
    case ( TypeIDs::int32_id ):
    {
      return lambda( int32(1) );
      break;
    }
    case ( TypeIDs::uint32_id ):
    {
      return lambda( uint32(1) );
      break;
    }
    case ( TypeIDs::int64_id ):
    {
      return lambda( int64(1) );
      break;
    }
    case ( TypeIDs::uint64_id ):
    {
      return lambda( uint64(1) );
      break;
    }
    case ( TypeIDs::real32_id ):
    {
      return lambda( real32(1) );
      break;
    }
    case ( TypeIDs::real64_id ):
    {
      return lambda( real64(1) );
      break;
    }
    case ( TypeIDs::r1Tensor_id ):
    {
      return lambda( R1Tensor() );
      break;
    }
    case ( TypeIDs::r2Tensor_id ):
    {
      return lambda( R2Tensor() );
      break;
    }
    case ( TypeIDs::r2SymTensor_id ):
    {
      return lambda( R2SymTensor() );
      break;
    }
    case ( TypeIDs::int32_array_id ):
    {
      return lambda( int32_array(1) );
      break;
    }
    case ( TypeIDs::uint32_array_id ):
    {
      return lambda( uint32_array(1) );
      break;
    }
    case ( TypeIDs::int64_array_id ):
    {
      return lambda( int64_array(1) );
      break;
    }
    case ( TypeIDs::uint64_array_id ):
    {
      return lambda( uint64_array(1) );
      break;
    }
    case ( TypeIDs::real32_array_id ):
    {
      return lambda( real32_array(1) );
      break;
    }
    case ( TypeIDs::real64_array_id ):
    {
      return lambda( real64_array(1) );
      break;
    }
    case ( TypeIDs::std_size_t_id ):
    {
      return lambda( std_size_t(1) );
      break;
    }
    case ( TypeIDs::string_id ):
    {
      return lambda( string("") );
      break;
    }
    default:
    {
      std::cout<<LOCATION<<std::endl;
      assert( false );
    }
    }
  }



  template< typename LAMBDA >
  static auto ApplyArrayTypeLambda1( const TypeIDs type,
                                         LAMBDA lambda )
  {
    switch( type )
    {
    case ( TypeIDs::int32_array_id ):
    {
      return lambda( int32_array(1) );
      break;
    }
    case ( TypeIDs::uint32_array_id ):
    {
      return lambda( uint32_array(1) );
      break;
    }
    case ( TypeIDs::int64_array_id ):
    {
      return lambda( int64_array(1) );
      break;
    }
    case ( TypeIDs::uint64_array_id ):
    {
      return lambda( uint64_array(1) );
      break;
    }
    case ( TypeIDs::real32_array_id ):
    {
      return lambda( real32_array(1) );
      break;
    }
    case ( TypeIDs::real64_array_id ):
    {
      return lambda( real64_array(1) );
      break;
    }
    case ( TypeIDs::r1_array_id ):
    {
      return lambda( r1_array(1) );
      break;
    }
    case ( TypeIDs::r2_array_id ):
    {
      return lambda( r2_array(1) );
      break;
    }
    case ( TypeIDs::r2Sym_array_id ):
    {
      return lambda( r2Sym_array(1) );
      break;
    }

    default:
    {
      std::cout<<LOCATION<<std::endl;
      assert( false );
    }
    }
  }

  template< typename LAMBDA >
  static auto ApplyTypeLambda1( const TypeIDs type,
                               LAMBDA lambda )
  {
    switch( type )
    {
    case ( TypeIDs::int32_id ):
    {
      return lambda( int32(1) );
      break;
    }
    case ( TypeIDs::uint32_id ):
    {
      return lambda( uint32(1) );
      break;
    }
    case ( TypeIDs::int64_id ):
    {
      return lambda( int64(1) );
      break;
    }
    case ( TypeIDs::uint64_id ):
    {
      return lambda( uint64(1) );
      break;
    }
    case ( TypeIDs::real32_id ):
    {
      return lambda( real32(1) );
      break;
    }
    case ( TypeIDs::real64_id ):
    {
      return lambda( real64(1) );
      break;
    }
    case ( TypeIDs::r1Tensor_id ):
    {
      return lambda( R1Tensor() );
      break;
    }
    case ( TypeIDs::r2Tensor_id ):
    {
      return lambda( R2Tensor() );
      break;
    }
    case ( TypeIDs::r2SymTensor_id ):
    {
      return lambda( R2SymTensor() );
      break;
    }
    case ( TypeIDs::int32_array_id ):
    {
      return lambda( int32_array(1) );
      break;
    }
    case ( TypeIDs::uint32_array_id ):
    {
      return lambda( uint32_array(1) );
      break;
    }
    case ( TypeIDs::int64_array_id ):
    {
      return lambda( int64_array(1) );
      break;
    }
    case ( TypeIDs::uint64_array_id ):
    {
      return lambda( uint64_array(1) );
      break;
    }
    case ( TypeIDs::real32_array_id ):
    {
      return lambda( real32_array(1) );
      break;
    }
    case ( TypeIDs::real64_array_id ):
    {
      return lambda( real64_array(1) );
      break;
    }
    case ( TypeIDs::r1_array_id ):
    {
      return lambda( r1_array(1) );
      break;
    }
    case ( TypeIDs::r2_array_id ):
    {
      return lambda( r2_array(1) );
      break;
    }
    case ( TypeIDs::r2Sym_array_id ):
    {
      return lambda( r2Sym_array(1) );
      break;
    }
    case ( TypeIDs::std_size_t_id ):
    {
      return lambda( std_size_t(1) );
      break;
    }
    case ( TypeIDs::string_id ):
    {
      return lambda( string("") );
      break;
    }
    case ( TypeIDs::string_array_id ):
    {
      return lambda( string_array(1) );
      break;
    }
    case ( TypeIDs::mapPair_array_id ):
    {
      return lambda( mapPair_array({1,1}) );
      break;
    }
    default:
    {
      std::cout<<LOCATION<<std::endl;
      assert( false );
    }
    }
  }


  template< typename LAMBDA >
  static auto ApplyTypeLambda2( const TypeIDs type,
                               LAMBDA lambda )
  {
    switch( type )
    {
    case ( TypeIDs::int32_id ):
    {
      return lambda( int32(1), int32(1) );
      break;
    }
    case ( TypeIDs::uint32_id ):
    {
      return lambda( uint32(1), uint32(1) );
      break;
    }
    case ( TypeIDs::int64_id ):
    {
      return lambda( int64(1), int64(1) );
      break;
    }
    case ( TypeIDs::uint64_id ):
    {
      return lambda( uint64(1), uint64(1) );
      break;
    }
    case ( TypeIDs::real32_id ):
    {
      return lambda( real32(1), real32(1) );
      break;
    }
    case ( TypeIDs::real64_id ):
    {
      return lambda( real64(1), real64(1) );
      break;
    }
    case ( TypeIDs::r1Tensor_id ):
    {
      return lambda( R1Tensor(), R1Tensor() );
      break;
    }
    case ( TypeIDs::r2Tensor_id ):
    {
      return lambda( R2Tensor(), R2Tensor() );
      break;
    }
    case ( TypeIDs::r2SymTensor_id ):
    {
      return lambda( R2SymTensor(), R2SymTensor() );
      break;
    }
    case ( TypeIDs::int32_array_id ):
    {
      return lambda( int32_array(1), int32(1) );
      break;
    }
    case ( TypeIDs::uint32_array_id ):
    {
      return lambda( uint32_array(1), uint32(1) );
      break;
    }
    case ( TypeIDs::int64_array_id ):
    {
      return lambda( int64_array(1), int64(1) );
      break;
    }
    case ( TypeIDs::uint64_array_id ):
    {
      return lambda( uint64_array(1), uint64(1) );
      break;
    }
    case ( TypeIDs::real32_array_id ):
    {
      return lambda( real32_array(1), real32(1) );
      break;
    }
    case ( TypeIDs::real64_array_id ):
    {
      return lambda( real64_array(1), real64(1) );
      break;
    }
    case ( TypeIDs::r1_array_id ):
    {
      return lambda( r1_array(1), R1Tensor() );
      break;
    }
    case ( TypeIDs::r2_array_id ):
    {
      return lambda( r2_array(1), R2Tensor() );
      break;
    }
    case ( TypeIDs::r2Sym_array_id ):
    {
      return lambda( r2Sym_array(1), R2SymTensor() );
      break;
    }
    case ( TypeIDs::std_size_t_id ):
    {
      return lambda( std_size_t(1), std_size_t(1) );
      break;
    }
    case ( TypeIDs::string_id ):
    {
      return lambda( string(""), string("") );
      break;
    }
    case ( TypeIDs::string_array_id ):
    {
      return lambda( string_array(1), string("") );
      break;
    }
//    case ( TypeIDs::mapPair_array_id ):
//    {
//      return lambda( mapPair_array(1), mapPair({}) );
//      break;
//    }

    default:
    {
      std::cout<<LOCATION<<std::endl;
      assert( false );
    }
    }
  }



  template< typename LAMBDA >
  static auto ApplyIntrinsicTypeLambda2( const TypeIDs type,
                                         LAMBDA lambda )
  {
    switch( type )
    {
    case ( TypeIDs::int32_id ):
    {
      return lambda( int32(1), int32(1) );
      break;
    }
    case ( TypeIDs::uint32_id ):
    {
      return lambda( uint32(1), uint32(1) );
      break;
    }
    case ( TypeIDs::int64_id ):
    {
      return lambda( int64(1), int64(1) );
      break;
    }
    case ( TypeIDs::uint64_id ):
    {
      return lambda( uint64(1), uint64(1) );
      break;
    }
    case ( TypeIDs::real32_id ):
    {
      return lambda( real32(1), real32(1) );
      break;
    }
    case ( TypeIDs::real64_id ):
    {
      return lambda( real64(1), real64(1) );
      break;
    }
    case ( TypeIDs::r1Tensor_id ):
    {
      return lambda( R1Tensor(), R1Tensor() );
      break;
    }
    case ( TypeIDs::r2Tensor_id ):
    {
      return lambda( R2Tensor(), R2Tensor() );
      break;
    }
    case ( TypeIDs::r2SymTensor_id ):
    {
      return lambda( R2SymTensor(), R2SymTensor() );
      break;
    }
    case ( TypeIDs::int32_array_id ):
    {
      return lambda( int32_array(1), int32(1) );
      break;
    }
    case ( TypeIDs::uint32_array_id ):
    {
      return lambda( uint32_array(1), uint32(1) );
      break;
    }
    case ( TypeIDs::int64_array_id ):
    {
      return lambda( int64_array(1), int64(1) );
      break;
    }
    case ( TypeIDs::uint64_array_id ):
    {
      return lambda( uint64_array(1), uint64(1) );
      break;
    }
    case ( TypeIDs::real32_array_id ):
    {
      return lambda( real32_array(1), real32(1) );
      break;
    }
    case ( TypeIDs::real64_array_id ):
    {
      return lambda( real64_array(1), real64(1) );
      break;
    }
    case ( TypeIDs::std_size_t_id ):
    {
      return lambda( std_size_t(1), std_size_t(1) );
      break;
    }
    case ( TypeIDs::string_id ):
    {
      return lambda( string(""), string("") );
      break;
    }
    case ( TypeIDs::string_array_id ):
    {
      return lambda( string_array(1), string("") );
      break;
    }
//    case ( TypeIDs::mapPair_array_id ):
//    {
//      return lambda( mapPair_array(1), mapPair({}) );
//      break;
//    }
    default:
    {
      std::cout<<LOCATION<<std::endl;
      assert( false );
    }
    }
  }







  inline static void equate( R1Tensor & lhs, int32 const component, real64 const & rhs )
  {
    lhs[component] = rhs;
  }

  template< typename TLHS, typename TRHS >
  inline static void equate( TLHS & lhs, int32 const component, TRHS const & rhs )
  {
    lhs = rhs;
  }


  inline static void add( R1Tensor & lhs, int32 const component, real64 const & rhs )
  {
    lhs[component] += rhs;
  }

  template< typename TLHS, typename TRHS >
  inline static void add( TLHS & lhs, int32 const component, TRHS const & rhs )
  {
    lhs += rhs;
  }



  inline static real64 value( R1Tensor & lhs, int32 const component )
  {
    return lhs[component];
  }

  inline static real64 value( R2Tensor & lhs, int32 const component )
  {
    return lhs.Data()[component];
  }

  inline static real64 value( R2SymTensor & lhs, int32 const component )
  {
    return lhs.Data()[component];
  }

  template< typename TLHS >
  inline static TLHS value( TLHS & lhs, int32 const component )
  {
    return lhs;
  }


  struct equateValue
  {
    inline static void f( R1Tensor & lhs, int32 const component, real64 const & rhs )
    {
      lhs[component] = rhs;
    }

    template< typename TLHS, typename TRHS >
    inline static void f( TLHS & lhs, int32 const component, TRHS const & rhs )
    {
      lhs = rhs;
    }

  };

  struct addValue
  {
    inline static void f( R1Tensor & lhs, int32 const component, real64 const & rhs )
    {
      lhs[component] += rhs;
    }

    inline static void f( R2Tensor & lhs, int32 const component, real64 const & rhs )
    {
      lhs.Data()[component] += rhs;
    }

    inline static void f( R2SymTensor & lhs, int32 const component, real64 const & rhs )
    {
      lhs.Data()[component] += rhs;
    }

    template< typename TLHS, typename TRHS >
    inline static void f( TLHS & lhs, int32 const component, TRHS const & rhs )
    {
      lhs += rhs;
    }

  };


  enum class operationType
  {
    add,
    multiply,
    equate
  };
};

}



#endif /* COMPONENTS_CORE_SRC_DATAREPOSITORY_DATATYPES_HPP_ */
