// Copyright (C) 2019 by Pedro Mendes, Rector and Visitors of the
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

// Copyright (C) 2005 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

/**
 * COptTask class.
 *
 * This class implements a optimization task which is comprised of a
 * of a problem and a method.
 *
 */
#include "copasi/copasi.h"

#include "COptTask.h"
#include "COptProblem.h"
#include "COptMethod.h"
#include "copasi/report/CKeyFactory.h"
#include "copasi/report/CReport.h"

#include "copasi/model/CModel.h"
#include "copasi/model/CState.h"

// static
const CTaskEnum::Method COptTask::ValidMethods[]  =
{
  CTaskEnum::Method::Statistics,
#ifdef COPASI_DEBUG
  CTaskEnum::Method::CoranaWalk,
#endif // COPASI_DEBUG
  CTaskEnum::Method::DifferentialEvolution,
  CTaskEnum::Method::SRES,
  CTaskEnum::Method::EvolutionaryProgram,
  CTaskEnum::Method::GeneticAlgorithm,
  CTaskEnum::Method::GeneticAlgorithmSR,
  CTaskEnum::Method::HookeJeeves,
  CTaskEnum::Method::LevenbergMarquardt,
  CTaskEnum::Method::NelderMead,
  CTaskEnum::Method::ParticleSwarm,
  CTaskEnum::Method::Praxis,
  CTaskEnum::Method::RandomSearch,
  CTaskEnum::Method::ScatterSearch,
  CTaskEnum::Method::SimulatedAnnealing,
  CTaskEnum::Method::SteepestDescent,
  CTaskEnum::Method::TruncatedNewton,
  CTaskEnum::Method::UnsetMethod
};

COptTask::COptTask(const CDataContainer * pParent,
                   const CTaskEnum::Task & type):
  CCopasiTask(pParent, type)
{
  mpProblem = new COptProblem(type, this);
  mpMethod = createMethod(CTaskEnum::Method::RandomSearch);

  ((COptMethod *) mpMethod)->setProblem((COptProblem *) mpProblem);
}

COptTask::COptTask(const COptTask & src,
                   const CDataContainer * pParent):
  CCopasiTask(src, pParent)
{
  mpProblem = new COptProblem(*(COptProblem *) src.mpProblem, this);
  mpMethod = createMethod(src.mpMethod->getSubType());
  //  mpMethod->setObjectParent(this);
  ((COptMethod *) mpMethod)->setProblem((COptProblem *) mpProblem);
}

COptTask::~COptTask()
{cleanup();}

void COptTask::cleanup() {}

bool COptTask::setCallBack(CProcessReport * pCallBack)
{
  bool success = CCopasiTask::setCallBack(pCallBack);

  if (!mpProblem->setCallBack(pCallBack)) success = false;

  if (!mpMethod->setCallBack(pCallBack)) success = false;

  return success;
}

bool COptTask::initialize(const OutputFlag & of,
                          COutputHandler * pOutputHandler,
                          std::ostream * pOstream)
{
  COptProblem * pProblem = dynamic_cast<COptProblem *>(mpProblem);
  COptMethod * pMethod = dynamic_cast<COptMethod *>(mpMethod);

  if (!pProblem || !pMethod) return false;

  //initialize reporting
  bool success = true;

  //do the part of the initialization of the subtask that needs to be
  //performed before the output is initialized. This is kind of a hack,
  //we need to find a more general solution for this
  success &= pProblem->initializeSubtaskBeforeOutput();

  success &= CCopasiTask::initialize(of, pOutputHandler, pOstream);

  //if (!mReport.open(pOstream)) success = false;
  //if (!mReport.compile()) success = false;

  success &= pProblem->initialize();

  pMethod->setProblem(pProblem);
  //  if (!pMethod->initialize()) return false;

  return success;
}

bool COptTask::process(const bool & useInitialValues)
{
  COptProblem * pProblem = dynamic_cast<COptProblem *>(mpProblem);
  COptMethod * pMethod = dynamic_cast<COptMethod *>(mpMethod);

  if (!pProblem || !pMethod) return false;

  mpMethod->isValidProblem(mpProblem);

  pProblem->randomizeStartValues();
  pProblem->rememberStartValues();

  if (useInitialValues) pProblem->resetEvaluations();

  output(COutputInterface::BEFORE);

  bool success = pMethod->optimise();

  pProblem->calculateStatistics();

  output(COutputInterface::AFTER);

  return success;
}

// virtual
const CTaskEnum::Method * COptTask::getValidMethods() const
{
  return COptTask::ValidMethods;
}
