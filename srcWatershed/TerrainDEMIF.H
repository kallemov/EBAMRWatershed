#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TERRAINDEMIF_H_
#define _TERRAINDEMIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
   This implicit function specifies a Digital Elevation Model
*/
class TerrainDEMIF: public BaseIF
{
public:
  ///
  /**
     Constructor specifying a digital elevation model
  */
    TerrainDEMIF(const ProblemDomain&     a_domain,
		 const int                a_interpType,
		 const RealVect&          a_finestDx,
		 const RealVect&          a_coarseDx,
		 const std::string&       a_demFile,
		 const Real               a_bottomElevation,
		 const RealVect           a_origin        = RealVect::Zero,
		 const Real               a_truncElev     = 1.e+123,
		 const Real               a_highGround    = 1.e+123,
		 const Real               a_verticalScale = 1.0,
		 const Real               a_cellvalue     = 1.0,
		 const int                a_edgeType      = 1);
    
    /// Copy constructor
    TerrainDEMIF(const TerrainDEMIF& a_inputIF);
    
    /// Destructor
    virtual ~TerrainDEMIF();

  ///
  /**
     Return the value of the function at a_point.
  */
  virtual Real value(const RealVect& a_point) const;
    
  virtual Real value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
		     const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

  /// Return true if DEM is read
  virtual bool readDEM(const std::string& a_demFile);

  ///
  /*
    Return true if DEM is read
  */
  virtual bool fixDEM();

protected:
  int      m_ncols;
  int      m_nrows;
  Real     m_xllcorner;
  Real     m_yllcorner;
  Real     m_cellvalue;
  Real     **m_DEM;
  RealVect m_finestDx;
  RealVect m_coarseDx;
  RealVect m_origin;
  IntVect  m_ncell;
  int      m_interpType;
  Real     m_highGround;
  Real     m_minDEM;
  Real     m_maxDEM;
  Real     m_bottomElevation;
  Real     m_truncElev;
  Real     m_verticalScale;
  bool     m_doCubic;
  Real     m_hx;
  Real     m_hy;
  Real     m_hx2;
  Real     m_hx3;
  Real     m_hy2;
  Real     m_hy3;
  int      m_edgeType;
private:
  TerrainDEMIF()
  {
    MayDay::Abort("TerrainDEMIF uses strong construction");
  }

  void operator=(const TerrainDEMIF& a_inputIF)
  {
    MayDay::Abort("TerrainDEMIF doesn't allow assignment");
  }

   Real getValueFromDEM(const int a_i, const int a_j) const;

};

#include "NamespaceFooter.H"
#endif
