// Copyright (C) 2019 - 2021 by Pedro Mendes, Rector and Visitors of the
// University of Virginia, University of Heidelberg, and University
// of Connecticut School of Medicine.
// All rights reserved.

// Copyright (C) 2017 - 2018 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and University of
// of Connecticut School of Medicine.
// All rights reserved.

// Copyright (C) 2010 - 2016 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 - 2009 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2002 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#ifdef WIN32
# pragma warning (disable: 4786)
# pragma warning (disable: 4243)
// warning C4355: 'this' : used in base member initializer list
# pragma warning (disable: 4355)
#endif  // WIN32

#include <limits.h>

#include <vector>
#include <numeric>
#include <limits>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <iostream>

#include "copasi/copasi.h"
#include "CHybridMethodRESD.h"
#include "copasi/core/CDataVector.h"
#include "copasi/function/CFunction.h"
#include "copasi/randomGenerator/CRandom.h"
#include "CTrajectoryMethod.h"
#include "CTrajectoryProblem.h"
#include "copasi/math/CMathContainer.h"
#include "copasi/model/CState.h"
#include "copasi/model/CCompartment.h"
#include "copasi/model/CModel.h"
#include "copasi/core/CRootContainer.h"
#include "copasi/model/CMetabNameInterface.h"


CHybridMethodRESD::CHybridMethodRESD(const CDataContainer * pParent,
                                       const CTaskEnum::Method & methodType,
                                       const CTaskEnum::Task & taskType):
  CTrajectoryMethod(pParent, methodType, taskType),
  mpRandomGenerator(NULL),
  mNumReactions(0),
  mMaxSteps(1000000),
  mNextReactionTime(0.0),
  mNextReactionIndex(C_INVALID_INDEX),
  mA0(0.0),
  mReactions(),
  mPropensityObjects(),
  mPropensityIdx(),
  mAmu(),
  mUpdateSequences(),
  mUpdateTimeDependentRoots(),
  mHaveTimeDependentRoots(false),
  mpRootValueCalculator(NULL),
  mMaxStepsReached(false),
  mTargetTime(0.0),
  mNumRoot(0),
  mRootsA(),
  mRootsB(),
  mRootsNonZero(),
  mpRootValueOld(NULL),
  mpRootValueNew(NULL),
  mapMetabIndex(),
  mThreshold(),
  switchAlgo(false),
  mDurationStochastic(0.0),
  mCountContainerVariables(0),
  rEmpiricalStandardDeviation(""),
  rEmpiricalStandardDeviationCalculated(false),
/*ODE45 variables*/
  mPropensitiesUpdateSequence(),
  mODE45(),
  mRKMethodStatus(CRungeKutta::INITIALIZE),
  mODEInitalized(false),
  mRootCounter(0),
  mMaxBalance(0),
  mData(),
  mY(),
  mpYdot(NULL),
  mSpeciesRateUpdateSequence(),
  mFirstReactionSpeciesIndex(C_INVALID_INDEX),
  mCountReactionSpecies(0),
  mContainerFluxes(),
  mpRelativeTolerance(NULL),
  mpAbsoluteTolerance(NULL),
  mpMaxInternalSteps(NULL),
  mLastRootTime(-std::numeric_limits< C_FLOAT64 >::infinity())

{
  assert((void *) &mData == (void *) &mData.dim);
  mData.pMethod = this;
  initializeParameter();
}

CHybridMethodRESD::CHybridMethodRESD(const CHybridMethodRESD & src,
                                       const CDataContainer * pParent):
  CTrajectoryMethod(src, pParent),
  mpRandomGenerator(NULL),
  mNumReactions(0),
  mMaxSteps(1000000),
  mNextReactionTime(0.0),
  mNextReactionIndex(C_INVALID_INDEX),
  mA0(0.0),
  mReactions(),
  mPropensityObjects(),
  mPropensityIdx(),
  mAmu(),
  mUpdateSequences(),
  mUpdateTimeDependentRoots(),
  mHaveTimeDependentRoots(false),
  mpRootValueCalculator(NULL),
  mMaxStepsReached(false),
  mTargetTime(src.mTargetTime),
  mNumRoot(src.mNumRoot),
  mRootsA(src.mRootsA),
  mRootsB(src.mRootsB),
  mRootsNonZero(src.mRootsNonZero),
  mpRootValueOld(NULL),
  mpRootValueNew(NULL),
  mapMetabIndex(),
  mThreshold(0.0),
  rEmpiricalStandardDeviation(""),
  switchAlgo(false),
  mCountContainerVariables(0),
  rEmpiricalStandardDeviationCalculated(false),
  mDurationStochastic(0.0)
{
  initializeParameter();
}

CHybridMethodRESD::~CHybridMethodRESD()
{
  if (mRootsFound.array() != NULL)
    {
      delete [] mRootsFound.array();
    }
}

void CHybridMethodRESD::initializeParameter()
{
  assertParameter("Duration for first stochastic", CCopasiParameter::Type::DOUBLE, (C_FLOAT64) 0.5);
  assertParameter("Relative Empirical Standard Deviation Threshold", CCopasiParameter::Type::DOUBLE, (C_FLOAT64) 1);
  /*-----------------------------------*/
  assertParameter("Max Internal Steps", CCopasiParameter::Type::INT, (C_INT32) 1000000);
  mpMaxInternalSteps = assertParameter("Max Internal Steps detr", CCopasiParameter::Type::UINT, (unsigned C_INT32) 100000);
  assertParameter("Use Random Seed", CCopasiParameter::Type::BOOL, (bool) false);
  assertParameter("Random Seed", CCopasiParameter::Type::UINT, (unsigned C_INT32) 1);
  //assertParameter("Runge Kutta Stepsize", CCopasiParameter::Type::DOUBLE, (C_FLOAT64) 0.1);
  mpRelativeTolerance = assertParameter("Relative Tolerance", CCopasiParameter::Type::UDOUBLE, (C_FLOAT64) 1.0e-006);
  mpAbsoluteTolerance = assertParameter("Absolute Tolerance", CCopasiParameter::Type::UDOUBLE, (C_FLOAT64) 1.0e-009);
  /* get configuration data */
  mMaxSteps = getValue< C_INT32 >("Max Internal Steps");
  /*get the Threshold*/
  mThreshold=getValue< C_FLOAT64 >("Relative Empirical Standard Deviation Threshold");
  /*get the first Stochastic Duration*/
  mDurationStochastic=getValue< C_FLOAT64 >("Duration for first stochastic");

  mpRootValueCalculator = new CBrent::EvalTemplate< CHybridMethodRESD >(this, &CHybridMethodRESD::rootValue);
}

bool CHybridMethodRESD::elevateChildren()
{
  initializeParameter();
  return true;
}

CTrajectoryMethod::Status CHybridMethodRESD::step(const double & deltaT,
    const bool & /* final */)
{
  C_FLOAT64 EndTime = *mpContainerStateTime + deltaT;

  if (mTargetTime != EndTime)
    {
      // We have a new end time and reset the root counter.
      mTargetTime = EndTime;
      mSteps = 0;
    }

  while (*mpContainerStateTime < EndTime)
    {
      // The Container State Time is updated during the reaction firing or root interpolation
      doSingleStep(*mpContainerStateTime, EndTime);

      if (mStatus == ROOT ||
          (mNumRoot > 0 && checkRoots()))
        {
          return ROOT;
        }

      if (mpProblem->getAutomaticStepSize())
        {
          break;
        }

      if (++mSteps > mMaxSteps)
        {
          CCopasiMessage(CCopasiMessage::EXCEPTION, MCTrajectoryMethod + 12);
        }
    }

  return NORMAL;
}

void CHybridMethodRESD::start()
{
  CTrajectoryMethod::start();

  //mStepsize = getValue< C_FLOAT64 >("Runge Kutta Stepsize");
  rEmpiricalStandardDeviationCalculated=false;
  mCountContainerVariables = mContainerState.size() - mpContainer->getCountFixedEventTargets();
  switchAlgo=false;
  /**
   * initialized mapMetabIndex
   * */
  size_t numModelValues=mpContainer->getModel().getModelValues().size();

  C_FLOAT64 numCopies=0.;
  for (size_t i = 0; i < numModelValues; ++i)
  {
   if(mpContainer->getModel().getModelValues()[i].getObjectName().compare("copy")==0)
   {
    numCopies=mpContainer->getModel().getModelValues()[i].getInitialValue();
    break;
   } 
  }
  numMetabs=mpContainer->getModel().getMetabolites().size();
  size_t numOriginalMetabs=numMetabs/(numCopies+1);

  mapMetabIndex.clear();
  mapMetabIndex.resize(numOriginalMetabs);
  C_FLOAT64 * pValue = mContainerState.begin();
  C_FLOAT64 * pValueEnd = mContainerState.end();
  CMathObject * pObject = mpContainer->getMathObject(pValue);
  rEmpiricalStandardDeviation+="Time\t\t";
  for (size_t i = 0; i < numOriginalMetabs; ++i)
  { 
   mapMetabIndex[i].first=mpContainer->getModel().getMetabolites()[i].getObjectName();
  rEmpiricalStandardDeviation+=mapMetabIndex[i].first+"\t\t";
  }
  rEmpiricalStandardDeviation+="\n";

  for (size_t i = 0; pValue != pValueEnd; ++pValue, ++pObject,++i)
  {
    std::string tmpNum=pObject->getObjectDisplayName().substr(0,pObject->getObjectDisplayName().find("."));

    for (size_t j = 0; j < numOriginalMetabs; ++j)
    {
      if (tmpNum.find("_")==0)
      {
        mapMetabIndex[j].second.push_back(i);
        break;
      }
      else 
        {
          if(mapMetabIndex[j].first.compare(tmpNum.substr(0,tmpNum.find("_")))==0)
          {
            mapMetabIndex[j].second.push_back(i);
            break;
          }
        }
    }
  }

  mFirstReactionSpeciesIndex = 1 + mpContainer->getCountODEs();
  mData.dim = mCountContainerVariables;
  mpRandomGenerator = &mpContainer->getRandomGenerator();
  mY.resize(mData.dim);
  mpYdot = mpContainer->getRate(false).array() + mpContainer->getCountFixedEventTargets();
  mCountReactionSpecies = mpContainer->getCountDependentSpecies() + mpContainer->getCountIndependentSpecies();
  mContainerFluxes.initialize(mpContainer->getFluxes());

 if (mRKMethodStatus==CRungeKutta::END)
  {
    mRKMethodStatus=CRungeKutta::RESTART;
  }
  CObjectInterface::ObjectSet Propensities;
  CObjectInterface::ObjectSet Fluxes;
   CMathReaction * itSlow = mReactions.begin();
   CMathReaction * endSlow = mReactions.end();
  for (; itSlow != endSlow; ++itSlow)
  {
    Propensities.insert(itSlow->getPropensityObject());
    Fluxes.insert(itSlow->getFluxObject());

  }
  mpContainer->getTransientDependencies().getUpdateSequence(mPropensitiesUpdateSequence, CCore::SimulationContext::Default, Fluxes,  Propensities);



  if (getValue< bool >("Use Random Seed"))
    {
      mpRandomGenerator->initialize(getValue< unsigned C_INT32 >("Random Seed"));
    }

  //mpCurrentState is initialized. This state is not used internally in the
  //stochastic solver, but it is used for returning the result after each step.

  //========Initialize Roots Related Arguments========
  mNumRoot = mpContainer->getRoots().size();

  if (mRootsFound.array() != NULL)
    {
      delete [] mRootsFound.array();
    }

  mRootsFound.initialize(mNumRoot, new C_INT[mNumRoot]);
  mRootsA.resize(mNumRoot);
  mRootsB.resize(mNumRoot);
  mpRootValueNew = &mRootsA;
  mpRootValueOld = &mRootsB;
  mRootsNonZero.resize(mNumRoot);
  mRootsNonZero = 0.0;
  mLastRootTime = -std::numeric_limits< C_FLOAT64 >::infinity();

  CMathObject * pRootObject = mpContainer->getMathObject(mpContainer->getRoots().array());
  CMathObject * pRootObjectEnd = pRootObject + mNumRoot;

  CObjectInterface::ObjectSet Requested;

  for (; pRootObject != pRootObjectEnd; ++pRootObject)
    {
      Requested.insert(pRootObject);
    }

  CObjectInterface::ObjectSet Changed;

  // Determine whether we have time dependent roots;

  CMathObject * pTimeObject = mpContainer->getMathObject(mpContainerStateTime);
  Changed.insert(pTimeObject);

  mpContainer->getTransientDependencies().getUpdateSequence(mUpdateTimeDependentRoots, CCore::SimulationContext::Default, Changed, Requested);
  mHaveTimeDependentRoots = (mUpdateTimeDependentRoots.size() > 0);

  // Build the reaction dependencies
  mReactions.initialize(mpContainer->getReactions());
  mNumReactions = mReactions.size();
  mAmu.initialize(mpContainer->getPropensities());
  mPropensityObjects.initialize(mAmu.size(), mpContainer->getMathObject(mAmu.array()));
  mUpdateSequences.resize(mNumReactions);

  CMathReaction * pReaction = mReactions.array();
  CMathReaction * pReactionEnd = pReaction + mNumReactions;
  CCore::CUpdateSequence * pUpdateSequence = mUpdateSequences.array();
  CMathObject * pPropensityObject = mPropensityObjects.array();
  CMathObject * pPropensityObjectEnd = pPropensityObject + mPropensityObjects.size();

  for (; pPropensityObject != pPropensityObjectEnd; ++pPropensityObject)
    {
      Requested.insert(pPropensityObject);
    }

  mPropensityIdx.resize(mNumReactions);
  size_t i = 0;

  for (; pReaction  != pReactionEnd; ++pReaction, ++pUpdateSequence, ++i)
    {
      Changed = pReaction->getChangedObjects();
      mPropensityIdx[i] = i;
      // The time is always updated
      Changed.insert(pTimeObject);

      pUpdateSequence->clear();
      mpContainer->getTransientDependencies().getUpdateSequence(*pUpdateSequence, CCore::SimulationContext::Default, Changed, Requested);
     
    }

  mMaxStepsReached = false;

  mTargetTime = *mpContainerStateTime;
  mNextReactionTime = *mpContainerStateTime;
  mNextReactionIndex = C_INVALID_INDEX;

  stateChange(CMath::eStateChange::State);
  
  return;
}

// virtual
bool CHybridMethodRESD::isValidProblem(const CCopasiProblem * pProblem)
{
  if (!CTrajectoryMethod::isValidProblem(pProblem)) return false;

  const CTrajectoryProblem * pTP = dynamic_cast<const CTrajectoryProblem *>(pProblem);

  if (pTP->getDuration() < 0.0)
    {
      //back integration not possible
      CCopasiMessage(CCopasiMessage::ERROR, MCTrajectoryMethod + 9);
      return false;
    }

  // check for ODEs
  if (mpContainer->getCountODEs() > 0)
    {
      CCopasiMessage(CCopasiMessage::ERROR, MCTrajectoryMethod + 28);
    }

  //TODO: rewrite CModel::suitableForStochasticSimulation() to use
  //      CCopasiMessage
  std::string message = mpContainer->getModel().suitableForStochasticSimulation();

  if (message != "")
    {
      //model not suitable, message describes the problem
      CCopasiMessage(CCopasiMessage::ERROR, message.c_str());
      return false;
    }

  if (getValue< C_INT32 >("Max Internal Steps") <= 0)
    {
      //max steps should be at least 1
      CCopasiMessage(CCopasiMessage::ERROR, MCTrajectoryMethod + 15);
      return false;
    }

  return true;
}

std::string CHybridMethodRESD::getrEmpiricalStandardDeviation() const
{
  return rEmpiricalStandardDeviation;
}

bool CHybridMethodRESD::calculteRelativeEmpiricalStandardDeviation(C_FLOAT64 startTime)
{
  std::string currentRelativeEmpiricalStandardDeviation;
  currentRelativeEmpiricalStandardDeviation+=std::to_string(startTime)+"\t";
  C_FLOAT64 * pValue = mContainerState.array() + mpContainer->getCountFixedEventTargets() + 1 /* Time */ + mpContainer->getCountODEs();
  /*we can get the number of Copies form initialized mapMetabIndex
  *size for vector postion
  */
  
  C_FLOAT64 average=0.;
  C_FLOAT64 variance=0.;
  C_FLOAT64 EmpiricalStandardDeviation=0.;
  size_t numOriginalMetabs=mapMetabIndex.size();
  size_t numCopies=mapMetabIndex[0].second.size();
  for (size_t i = 0; i < numOriginalMetabs; ++i)
  {
    
    average=0.;
    variance=0.;
    EmpiricalStandardDeviation=0.;
    for (size_t j = 0; j < numCopies; ++j)
      average+=pValue[mapMetabIndex[i].second[j]-1];
    
    average=average/numCopies;
    for (size_t j = 0; j < numCopies; ++j)
    {
      variance+=std::pow(pValue[mapMetabIndex[i].second[j]-1]-average,2.);
    }
    variance=variance/(numCopies-1);

    EmpiricalStandardDeviation=std::sqrt(variance)/average;
    currentRelativeEmpiricalStandardDeviation+=std::to_string(EmpiricalStandardDeviation)+"\t";
    vEmpiricalStandardDeviation.push_back(EmpiricalStandardDeviation);

  }
  currentRelativeEmpiricalStandardDeviation+="\n";
  rEmpiricalStandardDeviation+=currentRelativeEmpiricalStandardDeviation;
  rEmpiricalStandardDeviationCalculated=true;
#ifdef DEBUG_OUTPUT
  std::cout<<"max Threshold :"<<*max_element(vEmpiricalStandardDeviation.begin(), vEmpiricalStandardDeviation.end())<<"\n";
  std::cout<<"currentRelativeEmpiricalStandardDeviation :"<<currentRelativeEmpiricalStandardDeviation<<std::endl;
#endif // DEBUG_OUTPUT
  if (*max_element(vEmpiricalStandardDeviation.begin(), vEmpiricalStandardDeviation.end())<mThreshold)
  {
    switchAlgo=true;
    for (size_t i = 0; i < numOriginalMetabs; ++i)
    {
      
      average=0.;
      for (size_t j = 0; j < numCopies; ++j)
        average+=pValue[mapMetabIndex[i].second[j]-1];
      
      average=average/numCopies;
      for (size_t j = 0; j < numCopies; ++j)
        pValue[mapMetabIndex[i].second[j]-1]=average;
    }
  }
  return switchAlgo;
}

C_FLOAT64 CHybridMethodRESD::doSingleStep(C_FLOAT64 startTime, const C_FLOAT64 & endTime)
{
  if (!rEmpiricalStandardDeviationCalculated && startTime > mDurationStochastic )
      switchAlgo=calculteRelativeEmpiricalStandardDeviation(startTime);
  if (switchAlgo)
  {
#ifdef DEBUG_OUTPUT
    std::cout<<"old 0:"<<mContainerState<<"|"<<mContainerState.size()<<std::endl;
#endif // DEBUG_OUTPUT
    integrateDeterministicPart(endTime);
    *mpContainerStateTime = endTime;
    mStatus = NORMAL;
    return endTime - startTime;
  }
  else
  {

        if (mNextReactionIndex == C_INVALID_INDEX)
        {
          if (mA0 == 0)
            {
              *mpContainerStateTime = endTime;
              return endTime - startTime;
            }

          // We need to throw an exception if mA0 is NaN
          if (std::isnan(mA0))
            {
              CCopasiMessage(CCopasiMessage::EXCEPTION, MCTrajectoryMethod + 27);
            }

          mNextReactionTime = startTime - log(mpRandomGenerator->getRandomOO()) / mA0;

          // We are sure that we have at least 1 reaction
          C_FLOAT64 rand = mpRandomGenerator->getRandomOO() * mA0;
          size_t * idxProp = mPropensityIdx.begin();
          C_FLOAT64 sum = 0.0;
          size_t temp_prop;

          for (size_t i = 0; i != mNumReactions; ++idxProp, ++i)
            {
              sum += mAmu[*idxProp];
             
              if (sum > rand) break;

              if (i != 0 && mAmu[*idxProp] > mAmu[*(idxProp - 1)])
                std::swap(*idxProp, *(idxProp - 1));
            }
            

          if (idxProp == mPropensityIdx.end())
            --idxProp;

          mNextReactionIndex = *(idxProp);
        }

      *mpContainerStateTime = mNextReactionTime;

      // Check whether any time dependent root changes sign in curTime, mNextReactionTime.
      // If we need to interpolate return with that time value without firing any reaction

      if (mHaveTimeDependentRoots)
        {
          mpContainer->applyUpdateSequence(mUpdateTimeDependentRoots);

          if (checkRoots())
            {
              // Interpolate to find the first root
              C_FLOAT64 RootTime;
              C_FLOAT64 RootValue;
              CBrent::findRoot(startTime, mNextReactionTime, mpRootValueCalculator, &RootTime, &RootValue, 1e-9);

              if (RootTime > endTime)
                {
                  *mpContainerStateTime = endTime;
                  mpContainer->applyUpdateSequence(mUpdateTimeDependentRoots);
                  *mpRootValueNew = mpContainer->getRoots();
                  mStatus = NORMAL;

                  return endTime - startTime;
                }

              // If the last root time is equal to the current one we have already dealt with it and can proceed.
              // This is a save assumptions since reaction events advance the time and discrete events reset mLastRootTime.
              if (mLastRootTime < RootTime)
                {
                  mLastRootTime = RootTime;
                  *mpContainerStateTime = RootTime;
                  mpContainer->applyUpdateSequence(mUpdateTimeDependentRoots);
                  *mpRootValueNew = mpContainer->getRoots();

                  // Mark the appropriate root
                  C_INT * pRootFound = mRootsFound.array();
                  C_INT * pRootFoundEnd = pRootFound + mNumRoot;
                  C_FLOAT64 * pRootValue = mpRootValueNew->array();

                  for (; pRootFound != pRootFoundEnd; ++pRootFound, ++pRootValue)
                    if (*pRootValue == RootValue || *pRootValue == -RootValue)
                      {
                        *pRootFound = static_cast< C_INT >(CMath::RootToggleType::ToggleBoth);
                      }
                    else
                      {
                        *pRootFound = static_cast< C_INT >(CMath::RootToggleType::NoToggle);
                      }

                  mStatus = ROOT;

                  return RootTime - startTime;
                }
            }
        }

      if (mNextReactionTime >= endTime)
        {
          *mpContainerStateTime = endTime;
          mStatus = NORMAL;
          return endTime - startTime;
        }

      mReactions[mNextReactionIndex].fire();
      mpContainer->applyUpdateSequence(mUpdateSequences[mNextReactionIndex]);
      
      // calculate the total propensity
      mA0 = 0.0;

      const C_FLOAT64 * pAmu = mAmu.array();
      const C_FLOAT64 * pAmuEnd = pAmu + mNumReactions;

      for (; pAmu != pAmuEnd; ++pAmu)
        {
          mA0 += *pAmu;
        }

      mNextReactionIndex = C_INVALID_INDEX;
      mStatus = NORMAL;

      return mNextReactionTime - startTime; 
  }
}

/**
 * Check whether a root has been found
 */
bool CHybridMethodRESD::checkRoots()
{
  bool hasRoots = false;

  // Swap old an new for the next call.
  CVector< C_FLOAT64 > * pTmp = mpRootValueOld;
  mpRootValueOld = mpRootValueNew;
  mpRootValueNew = pTmp;

  *mpRootValueNew = mpContainer->getRoots();

  C_FLOAT64 *pRootValueOld = mpRootValueOld->array();
  C_FLOAT64 *pRootValueNew = mpRootValueNew->array();
  C_FLOAT64 *pRootNonZero = mRootsNonZero.array();

  C_INT *pRootFound    = mRootsFound.array();
  C_INT *pRootFoundEnd    = pRootFound + mRootsFound.size();

  const bool * pIsDiscrete = mpContainer->getRootIsDiscrete().array();
  const bool * pIsTimeDependent = mpContainer->getRootIsTimeDependent().array();

  for (; pRootFound != pRootFoundEnd; pRootValueOld++, pRootValueNew++, pRootFound++, pRootNonZero++, pIsDiscrete++, pIsTimeDependent++)
    {
      if (*pRootValueOld **pRootValueNew < 0.0 ||
          (*pRootValueNew == 0.0 && *pIsTimeDependent && !*pIsDiscrete))
        {
          // These root changes are not caused by the time alone as those are handled in do single step.
          hasRoots = true;
          *pRootFound = static_cast< C_INT >(CMath::RootToggleType::ToggleBoth);
        }
      else if (*pRootValueNew == 0.0 &&
               *pRootValueOld != 0.0)
        {
          hasRoots = true;
          *pRootFound = static_cast< C_INT >(CMath::RootToggleType::ToggleEquality); // toggle only equality
          *pRootNonZero = *pRootValueOld;
        }
      else if (*pRootValueNew != 0.0 &&
               *pRootValueOld == 0.0 &&
               *pRootValueNew **pRootNonZero < 0.0)
        {
          hasRoots = true;
          *pRootFound = static_cast< C_INT >(CMath::RootToggleType::ToggleInequality); // toggle only inequality
        }
      else
        {
          *pRootFound = static_cast< C_INT >(CMath::RootToggleType::NoToggle);
        }
    }

  return hasRoots;
}

/**
 * Update model state after one events happened
 */
void CHybridMethodRESD::stateChange(const CMath::StateChange & change)
{
  
  if (change & (CMath::StateChange(CMath::eStateChange::FixedEventTarget) | CMath::eStateChange::State | CMath::eStateChange::ContinuousSimulation | CMath::eStateChange::EventSimulation))
    {
      // Create a local copy of the state where the particle number species determined
      // by reactions are rounded to integers.

      C_FLOAT64 * pValue = mContainerState.array() + mpContainer->getCountFixedEventTargets() + 1 /* Time */ + mpContainer->getCountODEs();
      C_FLOAT64 * pValueEnd = pValue + mpContainer->getCountIndependentSpecies() + mpContainer->getCountDependentSpecies();

      for (; pValue != pValueEnd; ++pValue)
        {
          *pValue = floor(*pValue + 0.5);
        }

      // The container state is now up to date we just need to calculate all values needed for simulation.
      mpContainer->updateSimulatedValues(false); //for assignments

      CMathObject * pPropensityObject = mPropensityObjects.array();
      CMathObject * pPropensityObjectEnd = pPropensityObject + mPropensityObjects.size();
      C_FLOAT64 * pAmu = mAmu.array();
      mA0 = 0.0;

      // Update the propensity
      for (; pPropensityObject != pPropensityObjectEnd; ++pPropensityObject, ++pAmu)
        {
          pPropensityObject->calculateValue();
          mA0 += *pAmu;
        }

      mNextReactionIndex = C_INVALID_INDEX;
      *mpRootValueNew = mpContainer->getRoots();
      mLastRootTime = -std::numeric_limits< C_FLOAT64 >::infinity();
    }

  mMaxStepsReached = false;
  
  
}

C_FLOAT64 CHybridMethodRESD::rootValue(const C_FLOAT64 & time)
{
  *mpContainerStateTime = time;
  mpContainer->applyUpdateSequence(mUpdateTimeDependentRoots);

  const C_FLOAT64 * pRoot = mpContainer->getRoots().array();
  const C_FLOAT64 * pRootEnd = pRoot + mNumRoot;
  const C_FLOAT64 * pRootOld = mpRootValueOld->array();
  const C_FLOAT64 * pRootNew = mpRootValueNew->array();

  C_FLOAT64 MaxRootValue = - std::numeric_limits< C_FLOAT64 >::infinity();
  C_FLOAT64 RootValue;

  for (; pRoot != pRootEnd; ++pRoot, ++pRootOld, ++pRootNew)
    {
      // We are only looking for roots which change sign in [pOld, pNew]
      if (*pRootOld **pRootNew < 0 || *pRootNew == 0)
        {
          // Assure that the RootValue is increasing between old and new for each
          // candidate root.
          RootValue = (*pRootNew >= *pRootOld) ? *pRoot : -*pRoot;

          if (RootValue > MaxRootValue)
            {
              MaxRootValue = RootValue;
            }
        }
    }

  return MaxRootValue;
}




/**
 * Dummy f function for calculating derivative of y
 */
void CHybridMethodRESD::EvalF(const size_t * n, const C_FLOAT64 * t, const C_FLOAT64 * y,
                               C_FLOAT64 * ydot)
{static_cast<Data *>((void *) n)->pMethod->evalF(t, y, ydot);}

/**
 * Dummy f function for calculating roots value
 */
void CHybridMethodRESD::EvalR(const size_t * n, const C_FLOAT64 * t, const C_FLOAT64 * y,
                               const size_t * nr, C_FLOAT64 * r)
{static_cast<Data *>((void *) n)->pMethod->evalR(t, y, nr, r);}

/**
 * Derivative Calculation Function
 */
void CHybridMethodRESD::evalF(const C_FLOAT64 * t, const C_FLOAT64 * y, C_FLOAT64 * ydot)
{

  memcpy(mpContainerStateTime, y, mCountContainerVariables * sizeof(C_FLOAT64));
  *mpContainerStateTime = *t;

  mpContainer->updateSimulatedValues(false);
  memcpy(ydot, mpYdot, mCountContainerVariables * sizeof(C_FLOAT64));

  //(3) Deal with slow reactions
  CVector< C_FLOAT64 > mSavedFluxes;
      mSavedFluxes = mContainerFluxes;
      // Calculate the propensities of the slow reactions
      mpContainer->applyUpdateSequence(mPropensitiesUpdateSequence);

      // Calculate the species rates discarding the slow reaction fluxes.
      mpContainer->applyUpdateSequence(mSpeciesRateUpdateSequence);

      // Update species rates
      memcpy(ydot + mFirstReactionSpeciesIndex, mpYdot + mFirstReactionSpeciesIndex, mCountReactionSpecies * sizeof(C_FLOAT64));

      // Reset the fluxes
      mContainerFluxes = mSavedFluxes;
      mpContainer->applyUpdateSequence(mSpeciesRateUpdateSequence);
    

  return;
}

/**
 * Dummy Function for calculating roots value
 */
void CHybridMethodRESD::evalR(const C_FLOAT64 *t, const C_FLOAT64 *y,
                               const size_t *nr, C_FLOAT64 *r)
{

  memcpy(mpContainerStateTime, y, mCountContainerVariables * sizeof(C_FLOAT64));
  *mpContainerStateTime = *t;

  mpContainer->updateRootValues(false);

  CVectorCore< C_FLOAT64 > RootValues;


      RootValues.initialize(*nr - 1, r);

      C_FLOAT64 * pHybridRoot = r + (*nr - 1);

      const C_FLOAT64 * pAmu = y + mCountContainerVariables;
      const C_FLOAT64 * pAmuEnd = pAmu + mReactions.size();

      *pHybridRoot = mA0;

      for (; pAmu != pAmuEnd; ++pAmu)
        {
          *pHybridRoot -= *pAmu;
        }

  RootValues = mpContainer->getRoots();
  return;
}

//========Function for ODE45========
/**
 * Integrates the deterministic reactions of the system over
 * the specified time interval.
 *
 * @param ds A C_FLOAT64 specifying the stepsize.
 */
void CHybridMethodRESD::integrateDeterministicPart(C_FLOAT64 endTime)
{
  C_FLOAT64 StartTime = *mpContainerStateTime;

  //1----Set Parameters for ODE45 solver

  //=(3)= set time and old time

  //=(4)= set y and ode status
  if (mRKMethodStatus == CRungeKutta::INITIALIZE ||
      mRKMethodStatus == CRungeKutta::RESTART)
    {
      //only when starts a new step, we should copy state into ode solver
      memcpy(mY.array(), mpContainerStateTime, mCountContainerVariables * sizeof(C_FLOAT64));
    }
  else if (mRKMethodStatus != CRungeKutta::ERROR)
    {
      mRKMethodStatus = CRungeKutta::CONTINUE;
    }
  else
    {
      fatalError();
    }

  //3----If time increment is too small, do nothing
  C_FLOAT64 tdist = fabs(endTime - *mpContainerStateTime); //absolute time increment
  C_FLOAT64 w0 = std::max(fabs(*mpContainerStateTime), fabs(endTime));

  if (tdist < std::numeric_limits< C_FLOAT64 >::epsilon() * 2. * w0) //just do nothing
    {
      mRKMethodStatus = CRungeKutta::ERROR;
      *mpContainerStateTime = endTime;

      CCopasiMessage(CCopasiMessage::EXCEPTION, MCTrajectoryMethod + 6, "delta is too small");
      mRKMethodStatus = CRungeKutta::ERROR;

      return;
    }

  //4----just do nothing if there are no variables
  if (!mData.dim)
    {
      *mpContainerStateTime = endTime;
      return;
    }
    mY[0]=StartTime;
#ifdef DEBUG_OUTPUT
    std::cout<<"old 1:"<<mY<<"|"<<mY.size()<<std::endl;
#endif // DEBUG_OUTPUT

  //5----do integration
#ifdef DEBUG_OUTPUT
    std::cout<<" mpContainerStateTime :"<<*mpContainerStateTime<<std::endl;
#endif // DEBUG_OUTPUT
  mRKMethodStatus = mODE45(&mData.dim, mY.array(), mpContainerStateTime, &endTime,
                           mpContainer->getRoots().size(), (long int *)mpContainer->getRoots().array(), mRKMethodStatus, mpProblem->getAutomaticStepSize(),
                           mpRelativeTolerance, mpAbsoluteTolerance, mpMaxInternalSteps, EvalF, EvalR);
  
  C_FLOAT64 val=0.;
  
  size_t numOriginalMetabs=mapMetabIndex.size();
  size_t numCopies=mapMetabIndex[0].second.size();
  for (size_t i = 0; i < numOriginalMetabs; ++i)
  {
    
    val=mY[mapMetabIndex[i].second[numCopies-1]];
    for (size_t j = 0; j < numCopies; ++j)  
      mY[mapMetabIndex[i].second[j]]=val;
  }
  mY[0]=StartTime;
#ifdef DEBUG_OUTPUT
  std::cout<<"new 2:"<<mY<<"|"<<mY.size()<<std::endl;
#endif // DEBUG_OUTPUT
  
  // Update the container state and the simulated values.
    if (mRKMethodStatus!=CRungeKutta::ERROR)
    {
      memcpy(mpContainerStateTime, mY.array(), mCountContainerVariables * sizeof(C_FLOAT64));
      mpContainer->updateSimulatedValues(false);
    }

  return;
}