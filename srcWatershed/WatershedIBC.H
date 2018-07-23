#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _WATERSHEDIBC_H_
#define _WATERSHEDIBC_H_

#include "RichardsDomainBC.H"
#include "RichardsEBBC.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBISLayout.H"
#include "ProblemDomain.H"
#include "BaseEBBC.H"
#include "BaseDomainBC.H"
//#include "EBAMRWatershed.H"

#include "UsingNamespace.H"

class WatershedIBC
{
public:

    WatershedIBC(const bool  a_includeSurfaceSolver,
		 const IntVect&        a_scalarDomainLoBCType,
		 const IntVect&        a_scalarDomainHiBCType,
		 const Vector<RefCountedPtr<ParserBCValue> >& a_scalardomainLoBCValueFunc,
		 const Vector<RefCountedPtr<ParserBCValue> >& a_scalardomainHiBCValueFunc,
		 const RefCountedPtr<ParserBCValue>&   a_scalarEBBCValueFunc);

    virtual ~WatershedIBC();

    RefCountedPtr<RichardsDomainBCFactory> getScalarBC() const;
    RefCountedPtr<RichardsEBBCFactory>     getScalarEBBC() const;
    
    
    bool isSurfaceSolverOn()
    {
	return m_includeSurfaceSolver;
    }
    
    void setBoundarySource(Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > a_EBSource)
    {
	m_EBSource = a_EBSource;
    }

    Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& getBoundarySource()
    {
	return m_EBSource;
    }

protected:
    IntVect  m_scalarDomainLoBCType;
    IntVect  m_scalarDomainHiBCType;
    Vector<RefCountedPtr<ParserBCValue> >  m_scalarDomainLoBCValueFunc;
    Vector<RefCountedPtr<ParserBCValue> >  m_scalarDomainHiBCValueFunc;
    RefCountedPtr<ParserBCValue>   m_scalarEBBCValueFunc;

    int      m_scalarEBBCType;
    Real     m_scalarEBBCValue;
    //EBAMRWatershed* m_watershedPtr;
    bool m_includeSurfaceSolver;

private:
    Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > m_EBSource;
    //Vector<RefCountedPtr<LevelData<EBFluxFAB> > >  m_TempFlux;
    int m_finestLevel;
};

///
/**
 */
class WatershedIBCFactory
{
public:
    ///
    /**
     */
    WatershedIBCFactory(const IntVect&        a_scalarDomainLoBCType,
			const IntVect&        a_scalarDomainHiBCType,
			const Vector<RefCountedPtr<ParserBCValue> >& a_scalarDomainLoBCValueFunc,
			const Vector<RefCountedPtr<ParserBCValue> >& a_scalarDomainHiBCValueFunc,
			const RefCountedPtr<ParserBCValue>&   a_scalarEBBCValueFunc)
	{
	    m_scalarDomainLoBCType  = a_scalarDomainLoBCType;
	    m_scalarDomainHiBCType  = a_scalarDomainHiBCType;
	    m_scalarDomainLoBCValueFunc = a_scalarDomainLoBCValueFunc;
	    m_scalarDomainHiBCValueFunc = a_scalarDomainHiBCValueFunc;
	    m_scalarEBBCValueFunc       = a_scalarEBBCValueFunc;
	}

    ///
    virtual ~WatershedIBCFactory()
	{
	}

    virtual WatershedIBC * create(const bool  a_includeSurfaceSolver) const 
	{
	    WatershedIBC* ibc = new WatershedIBC(a_includeSurfaceSolver,
						 m_scalarDomainLoBCType, 
						 m_scalarDomainHiBCType, 
						 m_scalarDomainLoBCValueFunc, 
						 m_scalarDomainHiBCValueFunc,
						 m_scalarEBBCValueFunc);
    
	    return ibc;
	}

protected:
    IntVect  m_scalarDomainLoBCType;
    IntVect  m_scalarDomainHiBCType;
    Vector<RefCountedPtr<ParserBCValue> >  m_scalarDomainLoBCValueFunc;
    Vector<RefCountedPtr<ParserBCValue> >  m_scalarDomainHiBCValueFunc;
    RefCountedPtr<ParserBCValue>   m_scalarEBBCValueFunc;
};
#endif
