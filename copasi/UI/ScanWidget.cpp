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

// Copyright (C) 2003 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

//***  In this file I have put "//+++" in all places where something has to be added
//***  if a new scan item is introduced.

#include <QLabel>
#include <QPushButton>
#include <QLayout>
#include <QToolTip>
#include <QWhatsThis>
#include <QComboBox>
#include <QLineEdit>
#include <QGridLayout>
#include <QHBoxLayout>

#include <QInputDialog>

#include "copasi/copasi.h"

#include "ScanWidget.h"
#include "copasi/scan/CScanTask.h"
#include "copasi/scan/CScanProblem.h"
#include "copasi/scan/CScanMethod.h"
#include "copasi/CopasiDataModel/CDataModel.h"
#include "copasi/core/CRootContainer.h"

#include "CQTaskHeaderWidget.h"
#include "CQTaskBtnWidget.h"
#include "CQSimpleSelectionTree.h"
#include "CCopasiSelectionDialog.h"

#include "copasi/report/CKeyFactory.h"
#include "qtUtilities.h"
#include "CScanContainerWidget.h"
#include "copasi/utilities/CopasiTime.h"

//+++
//#include "CScanWidgetBreak.h"
#include "CScanWidgetRandom.h"
#include "CScanWidgetRepeat.h"
#include "CScanWidgetScan.h"
#include "CScanWidgetTask.h"

ScanWidget::ScanWidget(QWidget* parent, const char* name, Qt::WindowFlags f)
  : TaskWidget(parent, name, f)
{
  if (!name)
    setObjectName("ScanWidget");

  setWindowTitle("ScanWidget");
  ScanWidgetLayout = new QGridLayout(this);
#if QT_VERSION < QT_VERSION_CHECK(6,0,0)
  ScanWidgetLayout->setMargin(11);
#else
  ScanWidgetLayout->setContentsMargins(11, 11, 11, 11);
#endif
  ScanWidgetLayout->setSpacing(6);
  ScanWidgetLayout->setObjectName("ScanWidgetLayout");

  mpHeaderWidget->setTaskName("Parameter Scan");

  ScanWidgetLayout->addWidget(mpHeaderWidget, 0, 0);

  mpBtnWidget->verticalLayout->removeItem(mpBtnWidget->verticalSpacer);
  ScanWidgetLayout->addWidget(mpBtnWidget, 3, 0);

  //*****************

  QHBoxLayout* tmpLayout = new QHBoxLayout();

  QSpacerItem* tmpSpacer = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);
  tmpLayout->addItem(tmpSpacer);

  QLabel* tmpLabel = new QLabel(this);
  tmpLabel->setText("New scan item: ");
  tmpLayout->addWidget(tmpLabel);

  comboType = new QComboBox(this);
  //+++
  comboType->addItem("Scan");
  comboType->addItem("Repeat");
  comboType->addItem("Random distribution");
  tmpLayout->addWidget(comboType);

  QSpacerItem *mpSpacer = new QSpacerItem(20, 20, QSizePolicy::Maximum, QSizePolicy::Minimum);
  tmpLayout->addItem(mpSpacer);

  buttonNewItem = new QPushButton(this);
  buttonNewItem->setObjectName("buttonNewItem");
  buttonNewItem->setText("Create");
  //ScanWidgetLayout->addWidget(buttonNewItem, 1, 2);
  tmpLayout->addWidget(buttonNewItem);

  ScanWidgetLayout->addLayout(tmpLayout, 1, 0);

  //*****************************

  scrollview = new CScanContainerWidget(this);
  ScanWidgetLayout->addWidget(scrollview, 2, 0);

  // tab order
  /*setTabOrder(taskName, sExecutable);
  setTabOrder(sExecutable, steadyState);*/

  connect(buttonNewItem, SIGNAL(clicked()), this, SLOT(slotAddItem()));
}

ScanWidget::~ScanWidget()
{}

bool ScanWidget::runTask()
{
  if (!commonBeforeRunTask()) return false;

  bool success = true;

  if (!commonRunTask()) success = false;

  return success;
}

bool ScanWidget::loadTaskProtected()
{
  if (mIsLoading)
    return true; // currently loading

  mIsLoading = true;

  CScanTask* scanTask =
    dynamic_cast< CScanTask * >(mpObject);

  if (!scanTask) return false;

  CScanProblem *scanProblem = dynamic_cast<CScanProblem *>(scanTask->getProblem());

  if (!scanProblem) return false;

  scrollview->clearWidgetList();

  //std::vector<QWidget*> & widgetList = scrollview->getWidgetList();

  //+++
  CScanWidgetScan* tmp1;
  CScanWidgetRepeat* tmp2;
  CScanWidgetRandom* tmp3;
  //CScanWidgetBreak* tmp4;

  // the scan items
  assert(CRootContainer::getDatamodelList()->size() > 0);

  size_t i, imax = scanProblem->getNumberOfScanItems();

  for (i = 0; i < imax; ++i)
    {
      CScanProblem::Type type = scanProblem->getScanItem(i)->getValue< CScanProblem::Type >("Type");

      switch (type)
        {
          //+++
          case CScanProblem::SCAN_LINEAR:
            tmp1 = new CScanWidgetScan(scrollview);
            tmp1->load(scanProblem->getScanItem(i));
            scrollview->addWidget(tmp1);
            break;

          case CScanProblem::SCAN_REPEAT:
            tmp2 = new CScanWidgetRepeat(scrollview);
            tmp2->load(scanProblem->getScanItem(i));
            scrollview->addWidget(tmp2);
            break;

          case CScanProblem::SCAN_RANDOM:
            tmp3 = new CScanWidgetRandom(scrollview);
            tmp3->load(scanProblem->getScanItem(i));
            scrollview->addWidget(tmp3);
            break;

          default:
            break;
        }
    }

  // the widget for the subtask
  CScanWidgetTask* tmpT = new CScanWidgetTask(scrollview);
  tmpT->load(scanProblem);
  scrollview->addWidget(tmpT, false); //false: no control buttons (up/down/del)

  mChanged = false;

  mIsLoading = false;

  return true;
}

bool ScanWidget::slotAddItem()
{
  int totalRows = -1;
  //+++
  CScanWidgetScan* tmp1;
  CScanWidgetRepeat* tmp2;
  CScanWidgetRandom* tmp3;
  //CScanWidgetBreak* tmp4;

  int intType = comboType->currentIndex();
  CScanProblem::Type type;

  switch (intType)
    {
      case 0:
        type = CScanProblem::SCAN_LINEAR;
        break;

      case 1:
        type = CScanProblem::SCAN_REPEAT;
        break;

      case 2:
        type = CScanProblem::SCAN_RANDOM;
        break;

      default:
        type = CScanProblem::SCAN_LINEAR;
        break;
    }

  switch (type)
    {
      case CScanProblem::SCAN_REPEAT:
      {
        tmp2 = new CScanWidgetRepeat(scrollview);

        //create item to get the default values
        CCopasiParameterGroup* pItem = CScanProblem::createScanItem(type, 10);
        tmp2->load(pItem);
        pdelete(pItem);

        scrollview->insertWidget(tmp2);
        totalRows = scrollview->rowCount();
        tmp2->lineEditNumber->setFocus();
      }
      break;

      case CScanProblem::SCAN_LINEAR:
      {
        CQSimpleSelectionTree::ObjectClasses Classes = CQSimpleSelectionTree::InitialTime |
            CQSimpleSelectionTree::Parameters;

        std::vector< const CDataObject * > Selection = CCopasiSelectionDialog::getObjectVector(this, Classes);

        // create scan widgets as many as the number of selected objects
        std::vector< const CDataObject * >::iterator it = Selection.begin();
        std::vector< const CDataObject * >::iterator end = Selection.end();

        for (; it != end; ++it)
          {
            tmp1 = new CScanWidgetScan(scrollview);
            tmp1->initFromObject(*it);
            scrollview->insertWidget(tmp1);
            totalRows = scrollview->rowCount();
            tmp1->lineEditMin->setFocus();
          }

        break;
      }

      case CScanProblem::SCAN_RANDOM:
      {
        CQSimpleSelectionTree::ObjectClasses Classes = CQSimpleSelectionTree::InitialTime |
            CQSimpleSelectionTree::Parameters;

        std::vector< const CDataObject * > Selection = CCopasiSelectionDialog::getObjectVector(this, Classes);

        // create scan widgets as many as the number of selected objects
        std::vector< const CDataObject * >::iterator it = Selection.begin();
        std::vector< const CDataObject * >::iterator end = Selection.end();

        for (; it != end; ++it)
          {
            tmp3 = new CScanWidgetRandom(scrollview);
            tmp3->initFromObject(*it);
            scrollview->insertWidget(tmp3);
            totalRows = scrollview->rowCount();
            tmp3->lineEditMin->setFocus();
          }

        break;
      }

      default:
        break;
    }

  return true;
}

bool ScanWidget::saveTaskProtected()
{
  if (mIsLoading)
    return true; // currently loading the widget list cannot be trusted

  CScanTask* scanTask =
    dynamic_cast< CScanTask * >(mpObject);

  if (!scanTask) return false;

  CScanProblem *scanProblem = dynamic_cast<CScanProblem *>(scanTask->getProblem());

  if (!scanProblem) return false;

  const std::vector< QWidget * > & widgetList = scrollview->getWidgetList();

  size_t newSize = widgetList.size() - 1; // The last widget is reserved for the subtask.
  size_t oldSize = scanProblem->getNumberOfScanItems();

  size_t i, imax = std::min(newSize, oldSize);

  mChanged = false;

  for (i = 0; i < imax; ++i)
    {
      QWidget * pWidget = widgetList[i];

      if (pWidget->objectName() == "CScanWidgetScan")
        {
          mChanged |= static_cast< CScanWidgetScan * >(pWidget)->save(scanProblem->getScanItem(i));
        }
      else if (pWidget->objectName() == "CScanWidgetRandom")
        {
          mChanged |= static_cast< CScanWidgetRandom * >(pWidget)->save(scanProblem->getScanItem(i));
        }
      else if (pWidget->objectName() == "CScanWidgetRepeat")
        {
          mChanged |= static_cast< CScanWidgetRepeat * >(pWidget)->save(scanProblem->getScanItem(i));
        }
    }

  for (; i < newSize; ++i)
    {
      mChanged = true;
      QWidget * pWidget = widgetList[i];
      CCopasiParameterGroup * pItem = scanProblem->addScanItem(CScanProblem::SCAN_LINEAR);

      if (pWidget->objectName() == "CScanWidgetScan")
        {
          static_cast< CScanWidgetScan * >(pWidget)->save(pItem);
        }
      else if (pWidget->objectName() == "CScanWidgetRandom")
        {
          static_cast< CScanWidgetRandom * >(pWidget)->save(pItem);
        }
      else if (pWidget->objectName() == "CScanWidgetRepeat")
        {
          static_cast< CScanWidgetRepeat * >(pWidget)->save(pItem);
        }
    }

  for (; i < oldSize; ++i)
    {
      mChanged = true;
      scanProblem->removeScanItem(newSize);
    }

  // the subtask
  const CScanWidgetTask * tmpT = dynamic_cast<CScanWidgetTask*>(widgetList[newSize]);

  if (tmpT != NULL)
    {
      mChanged |= tmpT->save(scanProblem);
    }

  // :TODO Bug 322: This should only be called when actual changes have been saved.
  // However we do not check whether the scan item are mChanged we delete all
  // and create them new.
  if (mChanged)
    {
      if (mpDataModel != NULL)
        {
          mpDataModel->changed();
        }

      mChanged = false;
    }

  return true;
}

bool ScanWidget::taskFinishedEvent()
{
  if (!mpTask) return false;

  CScanProblem* pProblem = dynamic_cast<CScanProblem*>(mpTask->getProblem());

  if (!pProblem) return false;

  if (pProblem->getSubtask() != CTaskEnum::Task::parameterFitting)
    return false;

  protectedNotify(ListViews::ObjectType::MODELPARAMETERSET, ListViews::ADD);

  return true;
}

CCopasiMethod*ScanWidget::createMethod(const CTaskEnum::Method&) {return NULL;}
