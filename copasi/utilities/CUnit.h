// Copyright (C) 2019 - 2021 by Pedro Mendes, Rector and Visitors of the
// University of Virginia, University of Heidelberg, and University
// of Connecticut School of Medicine.
// All rights reserved.

// Copyright (C) 2017 - 2018 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and University of
// of Connecticut School of Medicine.
// All rights reserved.

// Copyright (C) 2014 - 2016 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// This class defines a unit, which might be considered
// an enhanced form of metadata, in COPASI. It can contains
// the expression of the unit, and a list of the SI units it is derived from.
// This is intended to facilitate unit analysis, user input parsing,
// and easy compatibility to the SBML units representation.

#ifndef CUNIT_H
#define CUNIT_H

#include "copasi/utilities/CUnitComponent.h"
#include "copasi/core/CDataVector.h"

class CUnitDefinition;

class CUnit
{
  friend std::ostream &operator<<(std::ostream &os, const CUnit & o);

public:

  typedef struct
  {
    std::string symbol;
    double multiplier;
    int scale;
    double exponent;
  } SymbolComponent;

  static C_FLOAT64 Avogadro;

  // Enum of unit types
  enum UnitType {time = 0, quantity, volume, area, length};

  // Enum of valid volume units
  enum VolumeUnit {dimensionlessVolume = 0, m3, l, ml, microl, nl, pl, fl};

  // String representation of valid volume units
  static const char * VolumeUnitNames[];

  // Enum of valid area units
  enum AreaUnit {dimensionlessArea = 0, m2, dm2, cm2, mm2, microm2, nm2, pm2, fm2};

  // String representation of valid area units
  static const char * AreaUnitNames[];

  // Enum of valid length units
  enum LengthUnit {dimensionlessLength = 0, m, dm, cm, mm, microm, nm, pm, fm};

  // String representation of valid length units
  static const char * LengthUnitNames[];

  // Enum of valid time units
  enum TimeUnit {dimensionlessTime = 0, d, h, min, s, ms, micros, ns, ps, fs, OldMinute};

  // String representation of valid time units
  static const char * TimeUnitNames[];

  //  Enum of valid quantity units
  enum QuantityUnit {dimensionlessQuantity = 0, Mol, mMol, microMol, nMol, pMol, fMol, number, OldXML};

  // String representation of valid quantity units as used in old (up to Build 18)
  // COPASI files
  static const char * QuantityUnitOldXMLNames[];

  // String representation of valid quantity units
  static const char * QuantityUnitNames[];

  static std::string replaceSymbol(const std::string & expression,
                                   const std::string & oldSymbol,
                                   const std::string & newSymbol);

  static std::string prettyPrint(const std::string & expression);

  // constructors
  /**
   * Default constructor
   */
  CUnit();

  /**
   * Kind constructor
   * @param const CBaseUnit::Kind & kind
   */
  CUnit(const CBaseUnit::Kind & kind);

  /**
   * Expression constructor
   * @param const std::string & expression
   */
  CUnit(const std::string & expression);

  /**
   * Copy constructor
   * @param const CUnit & src
   */
  CUnit(const CUnit & src);

  /**
   * Assignment operator
   */
  CUnit& operator=(const CUnit & rightSide);

  virtual ~CUnit();

  virtual CIssue setExpression(const std::string & expression);

  std::string getExpression() const;

  const std::set< std::string > & getUsedSymbols() const;

  void replaceSymbol(const std::string & oldSymbol,
                     const std::string & newSymbol);

  bool isDimensionless() const;
  bool isUndefined() const;
  bool isValid() const;

  void addComponent(const CUnitComponent & component);
  const std::set< CUnitComponent > & getComponents() const;

  CIssue compile();

  CUnit exponentiate(double exp) const;
  CUnit operator*(const CUnit & rightSide) const;
  bool operator==(const CUnit & rightSide) const;
  bool operator<(const CUnit & rightSide) const;
  bool isEquivalent(const CUnit & rightSide) const;
  void buildExpression();
  std::vector< SymbolComponent > getSymbolComponents() const;
  bool isUnitType(UnitType type) const;

  void setDimensionLess(double multiplier, double scale, double exponent);

private:
  std::string mExpression;
  std::set< CUnitComponent > mComponents;
  std::set< std::string > mUsedSymbols;
  CUnitComponent * mpDimensionless;

  static C_INT32 getExponentOfSymbol(const std::pair< std::string, CUnit > & SymbolDef, CUnit & unit);
  static std::pair< C_INT32, C_INT32 > removeSymbolFromUnit(const CUnit & symbol, CUnit & unit);
  static void replaceExponentInExpression(std::string & expression, const size_t & exponent);
};

#endif // CUNIT_H
