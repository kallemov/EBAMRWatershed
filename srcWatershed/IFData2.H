#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _IFDATA2_H_
#define _IFDATA2_H_

//This class is a copy of IFData class with modified constructor and initialization to support reuse
//this is done for effcient handling the DEM geometry case.


#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif


#include <map>
using std::map;
#include "Vector.H"
#include "REAL.H"
#include "IndexTM.H"
#include "IFData.H"
#include "Notation.H"
#include "IFSlicer.H"
#include "CoordinateSystem.H"

using std::endl;
using std::cerr;
#include "NamespaceHeader.H"

template <int dim> class IFData2
{
public:
  // data types to record vertex signs and edge intersections
  typedef IndexTM<int,dim>                      EdgeIndex;
  typedef map<EdgeIndex,Real > EdgeIntersections;

  typedef IndexTM<int,dim>                      Vertex;
  typedef map<Vertex,int >        CornerSigns;

  typedef IndexTM<int,dim>                      IvDim;
  typedef IndexTM<Real,dim>                     RvDim;

  typedef map<IvDim,RvDim >        NormalDerivatives;

  // member data
  CornerSigns           m_cornerSigns;
  EdgeIntersections     m_intersections;
  IFSlicer<dim>*        m_function;

  CoordinateSystem<dim> m_globalCoord;
  CoordinateSystem<dim> m_cellCenterCoord;
  CoordinateSystem<dim> m_parentCoord;
  CoordinateSystem<dim> m_localCoord;

  int                   m_maxOrder;
  NormalDerivatives     m_normalDerivatives;
  bool                  m_badNormal;

  bool                  m_allVerticesIn;
  bool                  m_allVerticesOut;
  bool                  m_allVerticesOn;
  bool                  m_isDefined;
  // empty constructor
  IFData2();

    // copy constructor
    IFData2(const IFData2<dim>& a_IFData);
    
    // makes an initial IFData2 from a BaseIF
    IFData2(const BaseIF & a_function);


  // makes an lower dimensional IFData2 from a IFSlicer; used for refinement (or it used to be)
    IFData2(IFSlicer<dim> * a_function);

  // destructor
  ~IFData2();

    void defineCoordinateSystem(const RvDim  & a_dx,
				const RvDim  & a_origin,
				const int    & a_maxOrder);

  void setNormalDerivatives();
  void makeCornerSigns();
  void findIntersectionPts();
  void defineLocalCoords();

  // used by the constructor
  bool isConnected(int          & a_edgeDir,
                   const Vertex & a_vertex1,
                   const Vertex & a_vertex2);

  void  makeEdgeKey(const int    & a_edgeDir,
                    const Vertex & a_vertex1,
                    const Vertex & a_vertex2);

  //this function calls the midPtRootFinder
  Real rootFinder(const EdgeIndex & a_thisEdge);

  Real BrentRootFinder(const RvDim & a_loPt,
                       const RvDim & a_hiPt,
                       const int   & a_edgeDir) const;
//                       const Real  & a_smallestRoot,
//                       const Real  & a_biggestRoot) const;

  void checkIntersection(bool       & a_hiOn,
                         bool       & a_loOn,
                         const Real & a_pt) const;

  void remakeCornerSigns();

  // output functions
  void print(ostream& out) const;

  // operators
  void operator=(const IFData2 & a_ifData);

  // this function sets the changing direction given the normal direction
  void setChangingDirection();

  bool searchNormalDir(const int i) const;
};

#include "NamespaceFooter.H"

#include "IFData2Implem.H"

#endif
