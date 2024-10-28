//
// Copyright (c) 2021, RTE (http://www.rte-france.com)
// See AUTHORS.txt
// All rights reserved.
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
//
// This file is part of Dynawo, an hybrid C++/Modelica open source time domain simulation tool for power systems.
//

/**
 * @file  DYNModelLoadOneTransformerTapChanger.h
 *
 * @brief ModelLoadOneTransformerTapChanger model header
 *
 */

#ifndef DYNAWO_SOURCES_MODELS_CPP_MODELLOADONETRANSFORMERTAPCHANGER_DYNMODELLOADONETRANSFORMERTAPCHANGER_H_
#define DYNAWO_SOURCES_MODELS_CPP_MODELLOADONETRANSFORMERTAPCHANGER_DYNMODELLOADONETRANSFORMERTAPCHANGER_H_

#include "DYNModelCPPImpl.h"
#include "DYNSubModelFactory.h"

namespace DYN {
class DataInterface;

/**
 * @brief ModelLoadOneTransformerTapChanger factory
 *
 * Implementation of @p SubModelFactory template for ModelLoadOneTransformerTapChanger Model
 */
class ModelLoadOneTransformerTapChangerFactory : public SubModelFactory {
 public:
  /**
   * @brief default constructor
   *
   */
  ModelLoadOneTransformerTapChangerFactory();

  /**
   * @brief default destructor
   *
   */
  ~ModelLoadOneTransformerTapChangerFactory();

  /**
   * @brief ModelLoadOneTransformerTapChanger getter
   *
   * @return A pointer to a new instance of ModelLoadOneTransformerTapChanger
   */
  SubModel* create() const;

  /**
   * @brief ModelLoadOneTransformerTapChanger destroy
   */
  void destroy(SubModel*) const;
};

/**
 * @brief ModelLoadOneTransformerTapChanger model class
 *
 *
 *
 */
class ModelLoadOneTransformerTapChanger : public ModelCPP::Impl {
 public:
  /**
   * @brief define type of calculated variables
   *
   */
  typedef enum {
    load_UPu = 0,
    transformer_P1Pu = 1,
    transformer_Q1Pu = 2,
    transformer_U1Pu = 3,
    tapChanger_valueMax = 4,
    tapChanger_valueMin = 5,
    state = 6,
    numCalculatedVars_ = 7
  } CalculatedVars_t;

  /**
   * @brief ModelLoadOneTransformerTapChanger model default constructor
   *
   *
   */
  ModelLoadOneTransformerTapChanger();

  /**
   * @brief ModelLoadOneTransformerTapChanger model default destructor
   *
   *
   */
  ~ModelLoadOneTransformerTapChanger();

  /**
   * @brief ModelLoadOneTransformerTapChanger model initialization
   * @param t0 : initial time of the simulation
   */
  void init(const double t0);

  /**
   * @brief ModelLoadOneTransformerTapChanger model's sizes getter
   *
   * Get the sizes of the vectors and matrices used by the solver to simulate
   * ModelLoadOneTransformerTapChanger instance. Used by @p ModelMulti to generate right size matrices
   * and vector for the solver.
   */
  void getSize();

  /**
   * @brief  ModelLoadOneTransformerTapChanger F(t,y,y') function evaluation
   *
   * Get the residues' values at a certain instant time with given state variables,
   * state variables derivatives
   * @param[in] t Simulation instant
   * @param[in] type type of the residues to compute (algebraic, differential or both)
   */
  void evalF(double t, propertyF_t type);

  /**
   * @brief ModelLoadOneTransformerTapChanger G(t,y,y') function evaluation
   *
   * Get the root's value
   *
   * @param t Simulation instant
   */
  void evalG(const double t);

  /**
   * @brief ModelLoadOneTransformerTapChanger discrete variables evaluation
   *
   * Get the discrete variables' value depending on current simulation instant and
   * current state variables values.
   *
   * @param t Simulation instant
   */
  void evalZ(const double t);

  /**
  * @brief set the silent flag for discrete variables
  * @param silentZTable flag table
  */
  void collectSilentZ(BitMask* silentZTable);

  /**
  * @brief Model mode change type evaluation
  *
  * Set the mode change type value depending on current simulation instant and
  * current state variables values.
  * @param[in] t Simulation instant
  * @return mode change type value
  */
  modeChangeType_t evalMode(const double t);

  /**
   * @brief ModelLoadOneTransformerTapChanger transposed jacobian evaluation
   *
   * Get the sparse transposed jacobian \f$ Jt=@F/@y + cj*@F/@y' \f$
   *
   * @param t Simulation instant
   * @param cj Jacobian prime coefficient
   * @param jt jacobian matrix to fullfill
   * @param rowOffset offset to use to identify the row where data should be added
   */
  void evalJt(const double t, const double cj, SparseMatrix& jt, const int rowOffset);

  /**
   * @brief  ModelLoadOneTransformerTapChanger transposed jacobian evaluation
   *
   * Get the sparse transposed jacobian \f$ Jt=@F/@y' \f$
   *
   * @param t Simulation instant
   * @param cj Jacobian prime coefficient
   * @param jt jacobian matrix to fullfill
   * @param rowOffset offset to use to identify the row where data should be added
   */
  void evalJtPrim(const double t, const double cj, SparseMatrix& jt, const int rowOffset);

  /**
   * @brief calculate calculated variables
   */
  void evalCalculatedVars();

  /**
   * @copydoc ModelCPP::getY0()
   */
  void getY0();

  /**
   * @copydoc ModelCPP::evalStaticYType()
   */
  void evalStaticYType();

  /**
   * @copydoc ModelCPP::evalDynamicYType()
   */
  void evalDynamicYType();

  /**
   * @copydoc ModelCPP::evalStaticFType()
   */
  void evalStaticFType();

  /**
   * @copydoc ModelCPP::evalDynamicFType()
   */
  void evalDynamicFType();

  /**
   * @brief get the index of variables used to define the jacobian associated to a calculated variable
   *
   * @param iCalculatedVar index of the calculated variable
   * @param indexes vector to fill with the indexes
   */
  void getIndexesOfVariablesUsedForCalculatedVarI(unsigned iCalculatedVar, std::vector<int>& indexes) const;

  /**
   * @brief evaluate the jacobian associated to a calculated variable
   *
   * @param iCalculatedVar index of the calculated variable
   * @param res values of the jacobian
   */
  void evalJCalculatedVarI(unsigned iCalculatedVar, std::vector<double>& res)const;

  /**
   * @brief evaluate the value of a calculated variable
   *
   * @param iCalculatedVar index of the calculated variable
   *
   * @return value of the calculated variable
   */
  double evalCalculatedVarI(unsigned iCalculatedVar) const;

  /**
   * @brief ModelLoadOneTransformerTapChanger parameters setter
   */
  void setSubModelParameters();

  /**
   * @brief ModelLoadOneTransformerTapChanger elements initializer
   *
   * Define elements for this model( elements to be seen by other models)
   *
   * @param elements  Reference to elements' vector
   * @param mapElement Map associating each element index in the elements vector to its name
   */

  void defineElements(std::vector<Element>& elements, std::map<std::string, int>& mapElement);

  /**
   * @brief initialize variables of the model
   *
   * A variable is a structure which contained all information needed to interact with the model
   * @param variables vector to fill with each variables
   */
  void defineVariables(std::vector<boost::shared_ptr<Variable> >& variables);

  /**
   * @brief define parameters
   * @param parameters vector to fill with each parameters
   */
  void defineParameters(std::vector<ParameterModeler>& parameters);

  /**
   * @brief get check sum number
   * @return the check sum number associated to the model
   */
  std::string getCheckSum() const;

  /**
   * @copydoc ModelCPP::initializeStaticData()
   */
  void initializeStaticData();

  /**
   * @brief initialize the model from data interface
   *
   * @param data data interface to use to initialize the model
   */
  void initializeFromData(const boost::shared_ptr<DataInterface>& data);

  /**
   * @copydoc ModelCPP::setFequations()
   */
  void setFequations();

  /**
   * @copydoc ModelCPP::setGequations()
   */
  void setGequations();

  /**
   * @copydoc ModelCPP::initParams()
   */
  void initParams();

  /**
  * @brief Coherence check on data (asserts, min/max values, sanity checks)
  */
  void checkDataCoherence(const double /*t*/) { /* not needed */ }

 private:
  static const int RUNNING_TRUE = 1;  ///< to represent running value
  static const int RUNNING_FALSE = 0;  ///< to represent a not running value

  static const int SWITCHOFF_TRUE = 1;  ///< to represent switchoff value
  static const int SWITCHOFF_FALSE = -1;  ///< to represent a not switchoff value

  static const int LOCKED_TRUE = 1;  ///< to represent locked value
  static const int LOCKED_FALSE = -1;  ///< to represent not locked value

  // Parameters
  double load_alpha_;  ///< Active load sensitivity to voltage
  double load_half_alpha_;  ///< 0.5 * alpha_
  double load_beta_;  ///< Reactive load sensitivity to voltage
  double load_half_beta_;  ///< 0.5 * beta_
  double transformer_B_;  ///< Susceptance in % (base U2Nom, SNom)
  double transformer_G_;  ///< Conductance in % (base U2Nom, SNom)
  double transformer_X_;  ///< Reactance in % (base U2Nom, SNom)
  double transformer_R_;  ///< Resistance in % (base U2Nom, SNom)
  int transformer_NbTap_;  ///< Number of taps
  double transformer_rTfoMaxPu_;  ///< Maximum transformation ratio in p.u: U2/U1 in no load conditions
  double transformer_rTfoMinPu_;  ///< Minimum transformation ratio in p.u: U2/U1 in no load conditions
  double transformer_P10Pu_;  ///< Start value of active power at terminal 1 in p.u (base SnRef) (receptor convention)
  double transformer_Q10Pu_;  ///< Start value of reactive power at terminal 1 in p.u (base SnRef) (receptor convention)
  double transformer_U10Pu_;  ///< Start value of voltage amplitude at terminal 1 in p.u (base UNom)
  double transformer_U1Phase0_;  ///< Start value of voltage angle at terminal 1 in rad
  double transformer_SNom_;  ///< Nominal apparent power in MVA

  double transformer_ZPu_re_;  ///< Transformer impedance in p.u (base U2Nom, SnRef) real part
  double transformer_ZPu_im_;  ///< Transformer impedance in p.u (base U2Nom, SnRef) imaginary part
  double transformer_YPu_re_;  ///< Transformer admittance in p.u (base U2Nom, SnRef) real part
  double transformer_YPu_im_;  ///< Transformer admittance in p.u (base U2Nom, SnRef) imaginary part

  double load_U0Pu_square_;  ///< Load's modulus of start value of voltage

  int running_;  ///< Indicates if the transformer is running or not at previous step
  double load_PRefPu_;  ///< Load's active power request at previous step
  double load_QRefPu_;  ///< Load's reactive power request at previous step
  int transformer_tap_;   ///< Transformer tap (between 0 and NbTap - 1) at previous step

  // Parameters
  double tapChanger_UDeadBand_;  ///< Voltage dead-band
  double tapChanger_UTarget_;  ///< Voltage set-point
  double tapChanger_t1st_;  ///< Time lag before changing the first tap
  double tapChanger_tNext_;  ///< Time lag before changing subsequent taps
  int tapChanger_tapMax_;  ///< Maximum tap
  int tapChanger_tapMin_;  ///< Minimum tap
  bool tapChanger_regulating0_;  ///< Whether the tap-changer is initially regulating

  // Internal parameters
  bool tapChanger_increaseTapToIncreaseValue_;  ///< Whether increasing the tap will increase the monitored value

  double whenUp_;  ///< when the voltage reached a value over the target+deadBand
  double whenDown_;     ///< when the voltage reached a value under the target-deadBand
  double whenLastTap_;  ///< last time when a tap changer
  bool moveUp_;         ///< @b true if tap should be increased
  bool moveDown_;       ///< @b false if tap should be decreased
  int tapRefDown_;      ///< initial tap when trying to decrease tap
  int tapRefUp_;        ///<  initial tap when trying to increase tap
  bool uMaxState_;      ///< @b true if U > uTarget + uDeadBand
  bool uMinState_;      ///< @b true if U < uTarget - uDeadBand
  bool uTargetState_;   ///< @b true if uTarget - uDeadBand < U < uTarget + uDeadBand

  /**
   * @brief enum to represent y indices value
   */
  typedef enum {
    transformer_terminal1_V_re = 0,
    transformer_terminal1_V_im = 1,
    load_V_re = 2,
    load_V_im = 3,
    load_i_re = 4,
    load_i_im = 5,
    load_PPu = 6,
    load_QPu = 7,
    transformer_terminal1_i_re = 8,
    transformer_terminal1_i_im = 9,
    numY = 10
  } yInd;

  /**
   * @brief enum to represent z indices value
   */
  typedef enum {
    load_PRefPu = 0,
    load_QRefPu = 1,
    transformer_rTfoPu = 2,
    transformer_tap = 3,
    tapChanger_locked = 4,
    running = 5,
    switchOffSignal1 = 6,
    switchOffSignal2 = 7,
    tapChanger_deltaUPuTarget = 8,
    numZ = 9
  } zInd;
};

}  // namespace DYN

#endif  // DYNAWO_SOURCES_MODELS_CPP_MODELLOADONETRANSFORMERTAPCHANGER_DYNMODELLOADONETRANSFORMERTAPCHANGER_H_
