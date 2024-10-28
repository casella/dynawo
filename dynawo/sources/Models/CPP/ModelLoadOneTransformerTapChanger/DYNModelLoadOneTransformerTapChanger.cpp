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
 * @file  DYNModelLoadOneTransformerTapChanger.cpp
 *
 * @brief ModelLoadOneTransformerTapChanger model implementation
 *
 */
#include <vector>

#include "DYNModelLoadOneTransformerTapChanger.h"
#include "DYNModelLoadOneTransformerTapChanger.hpp"
#include "DYNSparseMatrix.h"
#include "DYNElement.h"
#include "DYNCommonModeler.h"
#include "DYNVariableForModel.h"
#include "DYNModelConstants.h"
#include "DYNNumericalUtils.h"
// #include "DYNRTEMacrosMessage.h"

using std::vector;
using std::string;
using std::map;
using std::stringstream;

using boost::shared_ptr;

using parameters::ParametersSet;

extern "C" DYN::SubModelFactory* getFactory() {
  return (new DYN::ModelLoadOneTransformerTapChangerFactory());
}

extern "C" void deleteFactory(DYN::SubModelFactory* factory) {
  delete factory;
}

extern "C" DYN::SubModel* DYN::ModelLoadOneTransformerTapChangerFactory::create() const {
  DYN::SubModel* model(new DYN::ModelLoadOneTransformerTapChanger());
  return model;
}

extern "C" void DYN::ModelLoadOneTransformerTapChangerFactory::destroy(DYN::SubModel* model) const {
  delete model;
}

DYN::ModelLoadOneTransformerTapChangerFactory::ModelLoadOneTransformerTapChangerFactory() {
}

DYN::ModelLoadOneTransformerTapChangerFactory::~ModelLoadOneTransformerTapChangerFactory() {
}

namespace DYN {

ModelLoadOneTransformerTapChanger::ModelLoadOneTransformerTapChanger() :
Impl("ModelLoadOneTransformerTapChanger"),
load_alpha_(0.),
load_half_alpha_(0.),
load_beta_(0.),
load_half_beta_(0.),
transformer_B_(0.),
transformer_G_(0.),
transformer_X_(0.),
transformer_R_(0.),
transformer_NbTap_(0),
transformer_rTfoMaxPu_(0.),
transformer_rTfoMinPu_(0.),
transformer_P10Pu_(0.),
transformer_Q10Pu_(0.),
transformer_U10Pu_(0.),
transformer_U1Phase0_(0.),
transformer_SNom_(0.),
transformer_ZPu_re_(0.),
transformer_ZPu_im_(0.),
transformer_YPu_re_(0.),
transformer_YPu_im_(0.),
load_U0Pu_square_(0.),
running_(RUNNING_FALSE),
load_PRefPu_(0.),
load_QRefPu_(0.),
transformer_tap_(0),
tapChanger_UDeadBand_(0.),
tapChanger_UTarget_(0.),
tapChanger_t1st_(0.),
tapChanger_tNext_(0.),
tapChanger_tapMax_(0),
tapChanger_tapMin_(0),
tapChanger_regulating0_(false),
tapChanger_increaseTapToIncreaseValue_(false),
whenUp_(VALDEF),
whenDown_(VALDEF),
whenLastTap_(VALDEF),
moveUp_(false),
moveDown_(false),
tapRefDown_(-1),
tapRefUp_(-1),
uMaxState_(false),
uMinState_(false),
uTargetState_(true) {
}

ModelLoadOneTransformerTapChanger::~ModelLoadOneTransformerTapChanger() {
}

void
ModelLoadOneTransformerTapChanger::initializeFromData(const boost::shared_ptr<DataInterface>& /*data*/) {
}

void
ModelLoadOneTransformerTapChanger::initializeStaticData() {
}

void ModelLoadOneTransformerTapChanger::initParams() {
}

void
ModelLoadOneTransformerTapChanger::init(const double /*t0*/) {
}

void
ModelLoadOneTransformerTapChanger::getSize() {
  sizeF_ = numY - 2;  // 2 external variables
  sizeY_ = numY;
  sizeZ_ = numZ;
  sizeG_ = 6;
  sizeMode_ = 5;

  calculatedVars_.assign(numCalculatedVars_, 0);
}

void
ModelLoadOneTransformerTapChanger::evalStaticYType() {
  yType_[transformer_terminal1_V_re] = EXTERNAL;
  yType_[transformer_terminal1_V_im] = EXTERNAL;
  yType_[load_V_re] = ALGEBRAIC;
  yType_[load_V_im] = ALGEBRAIC;
  yType_[load_i_re] = ALGEBRAIC;
  yType_[load_i_im] = ALGEBRAIC;
  yType_[load_PPu] = ALGEBRAIC;
  yType_[load_QPu] = ALGEBRAIC;
  yType_[transformer_terminal1_i_re] = ALGEBRAIC;
  yType_[transformer_terminal1_i_im] = ALGEBRAIC;
}

void
ModelLoadOneTransformerTapChanger::evalStaticFType() {
  fType_[0] = ALGEBRAIC_EQ;
  fType_[1] = ALGEBRAIC_EQ;
  fType_[2] = ALGEBRAIC_EQ;
  fType_[3] = ALGEBRAIC_EQ;
  fType_[4] = ALGEBRAIC_EQ;
  fType_[5] = ALGEBRAIC_EQ;
  fType_[6] = ALGEBRAIC_EQ;
  fType_[7] = ALGEBRAIC_EQ;
}

void
ModelLoadOneTransformerTapChanger::evalDynamicYType() {
}

void
ModelLoadOneTransformerTapChanger::evalDynamicFType() {
}


void
ModelLoadOneTransformerTapChanger::evalF(double /*t*/, propertyF_t type) {
  if (type == UNDEFINED_EQ || type == ALGEBRAIC_EQ) {
    fLocal_[0] = yLocal_[load_PPu] - (yLocal_[load_V_re] * yLocal_[load_i_re] +
      yLocal_[load_V_im] * yLocal_[load_i_im]);
    fLocal_[1] = yLocal_[load_QPu] - (yLocal_[load_V_im] * yLocal_[load_i_re] -
      yLocal_[load_V_re] * yLocal_[load_i_im]);
    int runningVar = static_cast<int>(zLocal_[running]);
    if (runningVar == RUNNING_TRUE) {
      double load_UPu_over_U0Pu_square = (yLocal_[load_V_re] * yLocal_[load_V_re] +
        yLocal_[load_V_im] * yLocal_[load_V_im]) / load_U0Pu_square_;
      fLocal_[2] = yLocal_[load_PPu] - zLocal_[load_PRefPu] * pow_dynawo(load_UPu_over_U0Pu_square, load_half_alpha_);
      fLocal_[3] = yLocal_[load_QPu] - zLocal_[load_QRefPu] * pow_dynawo(load_UPu_over_U0Pu_square, load_half_beta_);
    } else {
      fLocal_[2] = yLocal_[load_i_re];
      fLocal_[3] = yLocal_[load_i_im];
    }

    if (runningVar == RUNNING_TRUE) {
      fLocal_[4] = yLocal_[transformer_terminal1_i_re] - zLocal_[transformer_rTfoPu] *
        (transformer_YPu_re_ * yLocal_[load_V_re] - transformer_YPu_im_ * yLocal_[load_V_im] + yLocal_[load_i_re]);
      fLocal_[5] = yLocal_[transformer_terminal1_i_im] - zLocal_[transformer_rTfoPu] * (transformer_YPu_re_ * yLocal_[load_V_im] +
        transformer_YPu_im_ * yLocal_[load_V_re] + yLocal_[load_i_im]);
      fLocal_[6] = transformer_ZPu_re_ * yLocal_[transformer_terminal1_i_re] - transformer_ZPu_im_ * yLocal_[transformer_terminal1_i_im] -
        (zLocal_[transformer_rTfoPu] * zLocal_[transformer_rTfoPu] * yLocal_[transformer_terminal1_V_re] -
        zLocal_[transformer_rTfoPu] * yLocal_[load_V_re]);
      fLocal_[7] = transformer_ZPu_im_ * yLocal_[transformer_terminal1_i_re] + transformer_ZPu_re_ * yLocal_[transformer_terminal1_i_im] -
        (zLocal_[transformer_rTfoPu] * zLocal_[transformer_rTfoPu] * yLocal_[transformer_terminal1_V_im] -
        zLocal_[transformer_rTfoPu] * yLocal_[load_V_im]);
    } else {
      fLocal_[4] = yLocal_[transformer_terminal1_i_re] - yLocal_[load_i_re];
      fLocal_[5] = yLocal_[transformer_terminal1_i_im] - yLocal_[load_i_im];
      fLocal_[6] = yLocal_[transformer_terminal1_V_re] - yLocal_[load_V_re];
      fLocal_[7] = yLocal_[transformer_terminal1_V_im] - yLocal_[load_V_im];
    }
  }
}

void
ModelLoadOneTransformerTapChanger::evalJt(const double /*t*/, const double /*cj*/, SparseMatrix& jt, const int rowOffset) {
  jt.changeCol();
  jt.addTerm(load_V_re + rowOffset, -yLocal_[load_i_re]);
  jt.addTerm(load_V_im + rowOffset, -yLocal_[load_i_im]);
  jt.addTerm(load_i_re + rowOffset, -yLocal_[load_V_re]);
  jt.addTerm(load_i_im + rowOffset, -yLocal_[load_V_im]);
  jt.addTerm(load_PPu + rowOffset, 1.);

  jt.changeCol();
  jt.addTerm(load_V_re + rowOffset, yLocal_[load_i_im]);
  jt.addTerm(load_V_im + rowOffset, -yLocal_[load_i_re]);
  jt.addTerm(load_i_re + rowOffset, -yLocal_[load_V_im]);
  jt.addTerm(load_i_im + rowOffset, yLocal_[load_V_re]);
  jt.addTerm(load_QPu + rowOffset, 1.);

  int runningVar = static_cast<int>(zLocal_[running]);
  if (runningVar == RUNNING_TRUE) {
    double load_UPu_square = yLocal_[load_V_re] * yLocal_[load_V_re] + yLocal_[load_V_im] * yLocal_[load_V_im];
    double load_UPu_over_U0Pu_square = load_UPu_square / load_U0Pu_square_;
    double derivative_coefficient1 = -zLocal_[load_PRefPu] * load_alpha_ * pow_dynawo(load_UPu_over_U0Pu_square, load_half_alpha_) /
      load_UPu_square;
    jt.changeCol();
    jt.addTerm(load_V_re + rowOffset, derivative_coefficient1 * yLocal_[load_V_re]);
    jt.addTerm(load_V_im + rowOffset, derivative_coefficient1 * yLocal_[load_V_im]);
    jt.addTerm(load_PPu + rowOffset, 1.);

    double derivative_coefficient2 = -zLocal_[load_QRefPu] * load_beta_ * pow_dynawo(load_UPu_over_U0Pu_square, load_half_beta_) /
      load_UPu_square;
    jt.changeCol();
    jt.addTerm(load_V_re + rowOffset, derivative_coefficient2 * yLocal_[load_V_re]);
    jt.addTerm(load_V_im + rowOffset, derivative_coefficient2 * yLocal_[load_V_im]);
    jt.addTerm(load_QPu + rowOffset, 1.);
  } else {
    jt.changeCol();
    jt.addTerm(load_i_re + rowOffset, 1.);

    jt.changeCol();
    jt.addTerm(load_i_im + rowOffset, 1.);
  }

  if (runningVar == RUNNING_TRUE) {
    jt.changeCol();
    jt.addTerm(load_V_re + rowOffset, -zLocal_[transformer_rTfoPu] * transformer_YPu_re_);
    jt.addTerm(load_V_im + rowOffset, zLocal_[transformer_rTfoPu] * transformer_YPu_im_);
    jt.addTerm(load_i_re + rowOffset, -zLocal_[transformer_rTfoPu]);
    jt.addTerm(transformer_terminal1_i_re + rowOffset, 1.);

    jt.changeCol();
    jt.addTerm(load_V_re + rowOffset, -zLocal_[transformer_rTfoPu] * transformer_YPu_im_);
    jt.addTerm(load_V_im + rowOffset, -zLocal_[transformer_rTfoPu] * transformer_YPu_re_);
    jt.addTerm(load_i_im + rowOffset, -zLocal_[transformer_rTfoPu]);
    jt.addTerm(transformer_terminal1_i_im + rowOffset, 1.);

    jt.changeCol();
    jt.addTerm(transformer_terminal1_V_re + rowOffset, -zLocal_[transformer_rTfoPu] * zLocal_[transformer_rTfoPu]);
    jt.addTerm(load_V_re + rowOffset, zLocal_[transformer_rTfoPu]);
    jt.addTerm(transformer_terminal1_i_re + rowOffset, transformer_ZPu_re_);
    jt.addTerm(transformer_terminal1_i_im + rowOffset, -transformer_ZPu_im_);

    jt.changeCol();
    jt.addTerm(transformer_terminal1_V_im + rowOffset, -zLocal_[transformer_rTfoPu] * zLocal_[transformer_rTfoPu]);
    jt.addTerm(load_V_im + rowOffset, zLocal_[transformer_rTfoPu]);
    jt.addTerm(transformer_terminal1_i_re + rowOffset, transformer_ZPu_im_);
    jt.addTerm(transformer_terminal1_i_im + rowOffset, transformer_ZPu_re_);

  } else {
    jt.changeCol();
    jt.addTerm(load_i_re + rowOffset, -1.);
    jt.addTerm(transformer_terminal1_i_re + rowOffset, 1.);

    jt.changeCol();
    jt.addTerm(load_i_im + rowOffset, -1.);
    jt.addTerm(transformer_terminal1_i_im + rowOffset, 1.);

    jt.changeCol();
    jt.addTerm(transformer_terminal1_V_re + rowOffset, 1.);
    jt.addTerm(load_V_re + rowOffset, -1.);

    jt.changeCol();
    jt.addTerm(transformer_terminal1_V_im + rowOffset, 1.);
    jt.addTerm(load_V_im + rowOffset, -1.);
  }
}

void
ModelLoadOneTransformerTapChanger::evalJtPrim(const double /*t*/, const double /*cj*/, SparseMatrix& jt, const int /*rowOffset*/) {
  jt.changeCol();
  jt.changeCol();
  jt.changeCol();
  jt.changeCol();
  jt.changeCol();
  jt.changeCol();
  jt.changeCol();
  jt.changeCol();
}

void
ModelLoadOneTransformerTapChanger::evalG(const double t) {
  double UMonitored = evalCalculatedVarI(load_UPu);
  int currentStepIndex = zLocal_[transformer_tap];
  int runningVar = static_cast<int>(zLocal_[running]);
  double lockedVar = zLocal_[tapChanger_locked];

  // ((switchoff_signal_1 = true or switchOffSignal2 = true) and running) || ((switchoff_signal_1 = false and switchOffSignal2 = false) && not(running))
  gLocal_[0] = (((zLocal_[switchOffSignal1] > 0 || zLocal_[switchOffSignal2] > 0) &&
    static_cast<int>(zLocal_[running]) == RUNNING_TRUE) || ((zLocal_[switchOffSignal1] < 0 && zLocal_[switchOffSignal2] < 0) &&
    static_cast<int>(zLocal_[running]) == RUNNING_FALSE)) ? ROOT_UP : ROOT_DOWN;

  gLocal_[1] = static_cast<int>(zLocal_[transformer_tap]) != transformer_tap_ ? ROOT_UP : ROOT_DOWN;

  double tapChanger_valueMax = tapChanger_UTarget_ + zLocal_[tapChanger_deltaUPuTarget] + tapChanger_UDeadBand_;
  double tapChanger_valueMin = tapChanger_UTarget_ + zLocal_[tapChanger_deltaUPuTarget] - tapChanger_UDeadBand_;
  gLocal_[2] = (UMonitored > tapChanger_valueMax && doubleNotEquals(UMonitored, tapChanger_valueMax) && runningVar == RUNNING_TRUE) ?
    ROOT_UP : ROOT_DOWN;  // U > Uc + deadBand
  gLocal_[3] = (UMonitored < tapChanger_valueMin && doubleNotEquals(UMonitored, tapChanger_valueMin) && runningVar == RUNNING_TRUE) ?
    ROOT_UP : ROOT_DOWN;  // U < Uc - deadBand
  gLocal_[4] = (moveUp_ && ((t - whenUp_ >= tapChanger_t1st_ && currentStepIndex == tapRefUp_) ||
    (t - whenLastTap_ >= tapChanger_tNext_ && currentStepIndex != tapRefUp_))
      && currentStepIndex < tapChanger_tapMax_ && !(lockedVar > 0) && runningVar == RUNNING_TRUE) ? ROOT_UP : ROOT_DOWN;  // first or next tap up
  gLocal_[5] = (moveDown_
    && ((t - whenDown_ >= tapChanger_t1st_ && currentStepIndex == tapRefDown_) || (t - whenLastTap_ >= tapChanger_tNext_ && currentStepIndex != tapRefDown_))
    && currentStepIndex > tapChanger_tapMin_ && !(lockedVar > 0) && runningVar == RUNNING_TRUE) ? ROOT_UP : ROOT_DOWN;  // first or next tap down
}

void
ModelLoadOneTransformerTapChanger::evalZ(const double t) {
  int runningVar = static_cast<int>(zLocal_[running]);
  int lockedVar = static_cast<int>(zLocal_[tapChanger_locked]);

  if (gLocal_[0] == ROOT_UP) {
    if (static_cast<int>(zLocal_[running]) == RUNNING_TRUE) {
      running_ = static_cast<int>(zLocal_[running]);
      zLocal_[running] = RUNNING_FALSE;
      DYNAddTimelineEvent(this, name(), LoadDisconnected);
    } else if (static_cast<int>(zLocal_[running]) == RUNNING_FALSE) {
      running_ = static_cast<int>(zLocal_[running]);
      zLocal_[running] = RUNNING_TRUE;
      DYNAddTimelineEvent(this, name(), LoadConnected);
    }
  }

  if (gLocal_[1] == ROOT_UP) {
    if (transformer_NbTap_ == 1)
      zLocal_[transformer_rTfoPu] = transformer_rTfoMinPu_;
    else
      zLocal_[transformer_rTfoPu] = transformer_rTfoMinPu_ + (transformer_rTfoMaxPu_ - transformer_rTfoMinPu_) * (zLocal_[transformer_tap] /
        (transformer_NbTap_ - 1));
  }

  if (runningVar == RUNNING_TRUE && !(lockedVar > 0)) {
    if (gLocal_[2] == ROOT_UP && !uMaxState_) {  // U > UMax
      if (!tapChanger_increaseTapToIncreaseValue_) {
        whenUp_ = t;
        moveUp_ = true;
        tapRefUp_ = static_cast<int>(zLocal_[transformer_tap]);
        whenDown_ = VALDEF;
        moveDown_ = false;
        tapRefDown_ = tapChanger_tapMax_;
      } else {
        whenDown_ = t;
        moveDown_ = true;
        tapRefDown_ = static_cast<int>(zLocal_[transformer_tap]);
        whenUp_ = VALDEF;
        moveUp_ = false;
        tapRefUp_ = tapChanger_tapMin_;
      }
      uMaxState_ = true;
      uMinState_ = false;
      uTargetState_ = false;
      DYNAddTimelineEvent(this, name(), TapChangerAboveMax);
    }

    if (gLocal_[3] == ROOT_UP && !uMinState_) {  // U< UMin
      if (!tapChanger_increaseTapToIncreaseValue_) {
        whenDown_ = t;
        moveDown_ = true;
        tapRefDown_ = static_cast<int>(zLocal_[transformer_tap]);
        whenUp_ = VALDEF;
        moveUp_ = false;
        tapRefUp_ = tapChanger_tapMin_;
      } else {
        whenUp_ = t;
        moveUp_ = true;
        tapRefUp_ = static_cast<int>(zLocal_[transformer_tap]);
        whenDown_ = VALDEF;
        moveDown_ = false;
        tapRefDown_ = tapChanger_tapMax_;
      }
      uMinState_ = true;
      uMaxState_ = false;
      uTargetState_ = false;
      DYNAddTimelineEvent(this, name(), TapChangerBelowMin);
    }

    if (gLocal_[2] == ROOT_DOWN && gLocal_[3] == ROOT_DOWN && !uTargetState_) {
      whenUp_ = VALDEF;
      moveUp_ = false;
      tapRefUp_ = tapChanger_tapMin_;
      whenDown_ = VALDEF;
      moveDown_ = false;
      tapRefDown_ = tapChanger_tapMax_;
      uTargetState_ = true;
      uMaxState_ = false;
      uMinState_ = false;
    }

    if (gLocal_[4] == ROOT_UP) {
      transformer_tap_ = static_cast<int>(zLocal_[transformer_tap]);
      zLocal_[transformer_tap] = zLocal_[transformer_tap] + 1;
      whenLastTap_ = t;
      DYNAddTimelineEvent(this, name(), TapUp);
    }

    if (gLocal_[5] == ROOT_UP) {
      transformer_tap_ = static_cast<int>(zLocal_[transformer_tap]);
      zLocal_[transformer_tap] = zLocal_[transformer_tap] - 1;
      whenLastTap_ = t;
      DYNAddTimelineEvent(this, name(), TapDown);
    }
  } else {
    whenUp_ = VALDEF;
    moveUp_ = false;
    tapRefUp_ = tapChanger_tapMin_;
    whenDown_ = VALDEF;
    moveDown_ = false;
    tapRefDown_ = tapChanger_tapMax_;
    uTargetState_ = false;
    uMaxState_ = false;
    uMinState_ = false;
    zLocal_[transformer_tap] = transformer_tap_;
  }
}

void
ModelLoadOneTransformerTapChanger::collectSilentZ(BitMask* silentZTable) {
  silentZTable[load_PRefPu].setFlags(NotUsedInDiscreteEquations);
  silentZTable[load_QRefPu].setFlags(NotUsedInDiscreteEquations);
  silentZTable[transformer_rTfoPu].setFlags(NotSilent);
  silentZTable[transformer_tap].setFlags(NotUsedInContinuousEquations);
  silentZTable[tapChanger_locked].setFlags(NotUsedInContinuousEquations);
  silentZTable[running].setFlags(NotSilent);
  silentZTable[switchOffSignal1].setFlags(NotUsedInContinuousEquations);
  silentZTable[switchOffSignal2].setFlags(NotUsedInContinuousEquations);
  silentZTable[tapChanger_deltaUPuTarget].setFlags(NotSilent);
}

modeChangeType_t
ModelLoadOneTransformerTapChanger::evalMode(const double /*t*/) {
  modeChangeType_t modeChangeType = NO_MODE;
  if (doubleNotEquals(load_PRefPu_, zLocal_[load_PRefPu])) {
    modeChangeType = ALGEBRAIC_MODE;
    load_PRefPu_ = zLocal_[load_PRefPu];
    // DYNRTEAddTimelineEvent(LoadTargetP, load_PRefPu_);
  }

  if (doubleNotEquals(load_QRefPu_, zLocal_[load_QRefPu])) {
    modeChangeType = ALGEBRAIC_MODE;
    load_QRefPu_ = zLocal_[load_QRefPu];
    // DYNRTEAddTimelineEvent(LoadTargetQ, load_QRefPu_);
  }

  if (transformer_tap_ != static_cast<int>(zLocal_[transformer_tap])) {
    modeChangeType = ALGEBRAIC_MODE;
    transformer_tap_ = static_cast<int>(zLocal_[transformer_tap]);
  }

  if (running_ != static_cast<int>(zLocal_[running])) {
    running_ = static_cast<int>(zLocal_[running]);
    return ALGEBRAIC_J_UPDATE_MODE;
  }

  return modeChangeType;
}

double
ModelLoadOneTransformerTapChanger::evalCalculatedVarI(unsigned iCalculatedVar) const {
  if (iCalculatedVar == load_UPu) {
    return std::sqrt(yLocal_[load_V_re] * yLocal_[load_V_re] + yLocal_[load_V_im] * yLocal_[load_V_im]);
  }

  if (iCalculatedVar == transformer_P1Pu) {
    return yLocal_[transformer_terminal1_V_re] * yLocal_[transformer_terminal1_i_re] +
      yLocal_[transformer_terminal1_V_im] * yLocal_[transformer_terminal1_i_im];
  }

  if (iCalculatedVar == transformer_Q1Pu) {
    return yLocal_[transformer_terminal1_V_im] * yLocal_[transformer_terminal1_i_re] -
      yLocal_[transformer_terminal1_V_re] * yLocal_[transformer_terminal1_i_im];
  }

  if (iCalculatedVar == transformer_U1Pu) {
    return std::sqrt(yLocal_[transformer_terminal1_V_re] * yLocal_[transformer_terminal1_V_re] +
      yLocal_[transformer_terminal1_V_im] * yLocal_[transformer_terminal1_V_im]);
  }

  if (iCalculatedVar == tapChanger_valueMax)
    return tapChanger_UTarget_ + zLocal_[tapChanger_deltaUPuTarget] + tapChanger_UDeadBand_;

  if (iCalculatedVar == tapChanger_valueMin)
    return tapChanger_UTarget_ + zLocal_[tapChanger_deltaUPuTarget] - tapChanger_UDeadBand_;

  if (iCalculatedVar == state) {
    if (static_cast<int>(zLocal_[running]) == RUNNING_TRUE) {
      return CLOSED;
    } else if (static_cast<int>(zLocal_[running]) == RUNNING_FALSE) {
      return OPEN;
    }
  }

  return 0.;
}

void
ModelLoadOneTransformerTapChanger::evalCalculatedVars() {
  calculatedVars_[load_UPu] = std::sqrt(yLocal_[load_V_re] * yLocal_[load_V_re] + yLocal_[load_V_im] * yLocal_[load_V_im]);

  calculatedVars_[transformer_P1Pu] = yLocal_[transformer_terminal1_V_re] * yLocal_[transformer_terminal1_i_re] +
    yLocal_[transformer_terminal1_V_im] * yLocal_[transformer_terminal1_i_im];
  calculatedVars_[transformer_Q1Pu] = yLocal_[transformer_terminal1_V_im] * yLocal_[transformer_terminal1_i_re] -
    yLocal_[transformer_terminal1_V_re] * yLocal_[transformer_terminal1_i_im];
  calculatedVars_[transformer_U1Pu] = std::sqrt(yLocal_[transformer_terminal1_V_re] * yLocal_[transformer_terminal1_V_re] +
    yLocal_[transformer_terminal1_V_im] * yLocal_[transformer_terminal1_V_im]);

  calculatedVars_[tapChanger_valueMax] = tapChanger_UTarget_ + zLocal_[tapChanger_deltaUPuTarget] + tapChanger_UDeadBand_;
  calculatedVars_[tapChanger_valueMin] = tapChanger_UTarget_ + zLocal_[tapChanger_deltaUPuTarget] - tapChanger_UDeadBand_;

  if (static_cast<int>(zLocal_[running]) == RUNNING_TRUE) {
    calculatedVars_[state] = CLOSED;
  } else if (static_cast<int>(zLocal_[running]) == RUNNING_FALSE) {
    calculatedVars_[state] = OPEN;
  }
}

void
ModelLoadOneTransformerTapChanger::getIndexesOfVariablesUsedForCalculatedVarI(unsigned iCalculatedVar, std::vector<int>& indexes) const {
  if (iCalculatedVar == load_UPu) {
    indexes.push_back(load_V_re);
    indexes.push_back(load_V_im);
  }

  if (iCalculatedVar == transformer_P1Pu) {
    indexes.push_back(transformer_terminal1_V_re);
    indexes.push_back(transformer_terminal1_V_im);
    indexes.push_back(transformer_terminal1_i_re);
    indexes.push_back(transformer_terminal1_i_im);
  }

  if (iCalculatedVar == transformer_Q1Pu) {
    indexes.push_back(transformer_terminal1_V_re);
    indexes.push_back(transformer_terminal1_V_im);
    indexes.push_back(transformer_terminal1_i_re);
    indexes.push_back(transformer_terminal1_i_im);
  }

  if (iCalculatedVar == transformer_U1Pu) {
    indexes.push_back(transformer_terminal1_V_re);
    indexes.push_back(transformer_terminal1_V_im);
  }
}

void
ModelLoadOneTransformerTapChanger::evalJCalculatedVarI(unsigned iCalculatedVar, vector<double>& res) const {
  if (iCalculatedVar == load_UPu) {
    double denominator = std::sqrt(yLocal_[load_V_re] * yLocal_[load_V_re] + yLocal_[load_V_im] * yLocal_[load_V_im]);
    res[0] = yLocal_[load_V_re] / denominator;
    res[1] = yLocal_[load_V_im] / denominator;
  }

  if (iCalculatedVar == transformer_P1Pu) {
    res[0] = yLocal_[transformer_terminal1_i_re];
    res[1] = yLocal_[transformer_terminal1_i_im];
    res[2] = yLocal_[transformer_terminal1_V_re];
    res[3] = yLocal_[transformer_terminal1_V_im];
  }

  if (iCalculatedVar == transformer_Q1Pu) {
    res[0] = -yLocal_[transformer_terminal1_i_im];
    res[1] = yLocal_[transformer_terminal1_i_re];
    res[2] = yLocal_[transformer_terminal1_V_im];
    res[3] = -yLocal_[transformer_terminal1_V_re];
  }

  if (iCalculatedVar == transformer_U1Pu) {
    double denominator = std::sqrt(yLocal_[transformer_terminal1_V_re] * yLocal_[transformer_terminal1_V_re] +
      yLocal_[transformer_terminal1_V_im] * yLocal_[transformer_terminal1_V_im]);
    res[0] = yLocal_[transformer_terminal1_V_re] / denominator;
    res[1] = yLocal_[transformer_terminal1_V_im] / denominator;
  }
}

void
ModelLoadOneTransformerTapChanger::getY0() {
  double transformer_puCoefficent = SNREF / (transformer_SNom_ * 100.);
  transformer_ZPu_re_ = transformer_R_ * transformer_puCoefficent;
  transformer_ZPu_im_ = transformer_X_ * transformer_puCoefficent;
  transformer_YPu_re_ = transformer_G_ * transformer_puCoefficent;
  transformer_YPu_im_ = transformer_B_ * transformer_puCoefficent;

  yLocal_[transformer_terminal1_V_re] = transformer_U10Pu_ * std::cos(transformer_U1Phase0_);
  yLocal_[transformer_terminal1_V_im] = transformer_U10Pu_ * std::sin(transformer_U1Phase0_);

  double transformer_UPu_square = yLocal_[transformer_terminal1_V_re] * yLocal_[transformer_terminal1_V_re] +
    yLocal_[transformer_terminal1_V_im] * yLocal_[transformer_terminal1_V_im];
  yLocal_[transformer_terminal1_i_re] = (transformer_P10Pu_ * yLocal_[transformer_terminal1_V_re] +
    transformer_Q10Pu_ * yLocal_[transformer_terminal1_V_im]) / transformer_UPu_square;
  yLocal_[transformer_terminal1_i_im] = (transformer_P10Pu_ * yLocal_[transformer_terminal1_V_im] -
    transformer_Q10Pu_ * yLocal_[transformer_terminal1_V_re]) / transformer_UPu_square;

  double deltauPu_re = transformer_ZPu_re_ * yLocal_[transformer_terminal1_i_re] - transformer_ZPu_im_ * yLocal_[transformer_terminal1_i_im];
  double deltauPu_im = transformer_ZPu_im_ * yLocal_[transformer_terminal1_i_re] + transformer_ZPu_re_ * yLocal_[transformer_terminal1_i_im];

  if (transformer_NbTap_ == 1) {
    zLocal_[transformer_tap] = 0;
    zLocal_[transformer_rTfoPu] = transformer_rTfoMinPu_;
  } else {
    double transformer_u10Pu = std::sqrt(yLocal_[transformer_terminal1_V_re] * yLocal_[transformer_terminal1_V_re] +
      yLocal_[transformer_terminal1_V_im] * yLocal_[transformer_terminal1_V_im]);
    if (doubleIsZero(transformer_u10Pu)) {
      zLocal_[transformer_tap] = 0;
      zLocal_[transformer_rTfoPu] = transformer_rTfoMinPu_;
    } else {
      double a = transformer_u10Pu * transformer_u10Pu;
      double b = -2. * deltauPu_re * yLocal_[transformer_terminal1_V_re] - 2. * deltauPu_im * yLocal_[transformer_terminal1_V_im] -
        tapChanger_UTarget_ * tapChanger_UTarget_;
      double deltauPuModule = std::sqrt(deltauPu_re * deltauPu_re + deltauPu_im * deltauPu_im);
      double c = deltauPuModule * deltauPuModule;
      double delta = b * b - 4. * a * c;

//      if (delta < 0.)
//        DYNRTEError(Error::MODELER, InvalidTransformerInitDelta, name(), delta);
      double root = (-b + sqrt(delta)) / (2. * a);
//      if (root < 0.)
//        DYNError(Error::MODELER, InvalidTransformerInitRoot, name(), root);
      double rcTfo0Pu = sqrt(root);

      double tapEstimation = ((rcTfo0Pu - transformer_rTfoMinPu_) / (transformer_rTfoMaxPu_ - transformer_rTfoMinPu_)) * (transformer_NbTap_ - 1);

      int tap0;
      if (tapEstimation < 0. || doubleIsZero(tapEstimation))
        tap0 = 0;
      else if (tapEstimation > (transformer_NbTap_ - 1) || doubleEquals(tapEstimation, transformer_NbTap_ - 1))
        tap0 = transformer_NbTap_ - 1;
      else if ((tapEstimation - std::floor(tapEstimation)) < (std::ceil(tapEstimation) - tapEstimation))
        tap0 = static_cast<int>(std::floor(tapEstimation));
      else
        tap0 = static_cast<int>(std::ceil(tapEstimation));
      zLocal_[transformer_tap] = tap0;
      zLocal_[transformer_rTfoPu] = transformer_rTfoMinPu_ + (transformer_rTfoMaxPu_ - transformer_rTfoMinPu_) * (zLocal_[transformer_tap] /
        (transformer_NbTap_ - 1));
    }
  }

  transformer_tap_ = zLocal_[transformer_tap];

  double u20Pu_re = (zLocal_[transformer_rTfoPu] * zLocal_[transformer_rTfoPu] * yLocal_[transformer_terminal1_V_re] - deltauPu_re) /
    zLocal_[transformer_rTfoPu];
  double u20Pu_im = (zLocal_[transformer_rTfoPu] * zLocal_[transformer_rTfoPu] * yLocal_[transformer_terminal1_V_im] - deltauPu_im) /
    zLocal_[transformer_rTfoPu];

  double YPu_u20Pu_re = transformer_YPu_re_ * u20Pu_re - transformer_YPu_im_ * u20Pu_im;
  double YPu_u20Pu_im = transformer_YPu_im_ * u20Pu_re + transformer_YPu_re_ * u20Pu_im;

  double i20Pu_re = zLocal_[transformer_rTfoPu] * YPu_u20Pu_re - yLocal_[transformer_terminal1_i_re] / zLocal_[transformer_rTfoPu];
  double i20Pu_im = zLocal_[transformer_rTfoPu] * YPu_u20Pu_im - yLocal_[transformer_terminal1_i_im] / zLocal_[transformer_rTfoPu];

  load_U0Pu_square_ = u20Pu_re * u20Pu_re + u20Pu_im * u20Pu_im;
  yLocal_[load_V_re] = u20Pu_re;
  yLocal_[load_V_im] = u20Pu_im;
  yLocal_[load_i_re] = -i20Pu_re;
  yLocal_[load_i_im] = -i20Pu_im;
  yLocal_[load_PPu] = yLocal_[load_V_re] * yLocal_[load_i_re] + yLocal_[load_V_im] * yLocal_[load_i_im];
  yLocal_[load_QPu] = yLocal_[load_V_im] * yLocal_[load_i_re] - yLocal_[load_V_re] * yLocal_[load_i_im];

  running_ = RUNNING_TRUE;
  load_PRefPu_ = yLocal_[load_PPu];
  load_QRefPu_ = yLocal_[load_QPu];
  zLocal_[load_PRefPu] = yLocal_[load_PPu];
  zLocal_[load_QRefPu] = yLocal_[load_QPu];
  zLocal_[running] = running_;

  zLocal_[switchOffSignal1] = SWITCHOFF_FALSE;
  zLocal_[switchOffSignal2] = SWITCHOFF_FALSE;

  zLocal_[tapChanger_deltaUPuTarget] = 0.;

  zLocal_[tapChanger_locked] = LOCKED_FALSE;
}

void
ModelLoadOneTransformerTapChanger::defineVariables(vector<shared_ptr<Variable> >& variables) {
  variables.push_back(VariableNativeFactory::createState("transformer_terminal1_V_re", CONTINUOUS));
  variables.push_back(VariableNativeFactory::createState("transformer_terminal1_V_im", CONTINUOUS));

  variables.push_back(VariableNativeFactory::createState("load_V_re", CONTINUOUS));
  variables.push_back(VariableNativeFactory::createState("load_V_im", CONTINUOUS));
  variables.push_back(VariableNativeFactory::createState("load_i_re", FLOW));
  variables.push_back(VariableNativeFactory::createState("load_i_im", FLOW));
  variables.push_back(VariableNativeFactory::createState("load_PPu", CONTINUOUS));
  variables.push_back(VariableNativeFactory::createState("load_QPu", CONTINUOUS));

  variables.push_back(VariableNativeFactory::createState("transformer_terminal1_i_re", FLOW));
  variables.push_back(VariableNativeFactory::createState("transformer_terminal1_i_im", FLOW));

  variables.push_back(VariableNativeFactory::createState("load_PRefPu", DISCRETE));
  variables.push_back(VariableNativeFactory::createState("load_QRefPu", DISCRETE));

  variables.push_back(VariableNativeFactory::createCalculated("load_UPu", CONTINUOUS));

  variables.push_back(VariableAliasFactory::create("transformer_terminal2_V_re", "load_V_re", CONTINUOUS));
  variables.push_back(VariableAliasFactory::create("transformer_terminal2_V_im", "load_V_im", CONTINUOUS));
  variables.push_back(VariableAliasFactory::create("transformer_terminal2_i_re", "load_i_re", FLOW, true));
  variables.push_back(VariableAliasFactory::create("transformer_terminal2_i_im", "load_i_im", FLOW, true));

  variables.push_back(VariableNativeFactory::createState("transformer_rTfoPu", DISCRETE));
  variables.push_back(VariableNativeFactory::createState("transformer_tap", DISCRETE));

  variables.push_back(VariableNativeFactory::createCalculated("transformer_P1Pu", CONTINUOUS));
  variables.push_back(VariableNativeFactory::createCalculated("transformer_Q1Pu", CONTINUOUS));
  variables.push_back(VariableNativeFactory::createCalculated("transformer_U1Pu", CONTINUOUS));
  variables.push_back(VariableAliasFactory::create("transformer_U2Pu", "load_UPu", CONTINUOUS));

  variables.push_back(VariableNativeFactory::createState("tapChanger_locked", BOOLEAN));

  variables.push_back(VariableNativeFactory::createState("running", BOOLEAN));
  variables.push_back(VariableNativeFactory::createState("switchOffSignal1", BOOLEAN));
  variables.push_back(VariableNativeFactory::createState("switchOffSignal2", BOOLEAN));

  variables.push_back(VariableNativeFactory::createState("tapChanger_deltaUPuTarget", DISCRETE));

  variables.push_back(VariableNativeFactory::createCalculated("tapChanger_valueMax", CONTINUOUS));
  variables.push_back(VariableNativeFactory::createCalculated("tapChanger_valueMin", CONTINUOUS));

  variables.push_back(VariableNativeFactory::createCalculated("state", INTEGER));
}

void
ModelLoadOneTransformerTapChanger::defineElements(std::vector<Element>& elements, std::map<std::string, int>& mapElement) {
  addElement("transformer", Element::STRUCTURE, elements, mapElement);
  addSubElement("terminal1", "transformer", Element::STRUCTURE, name(), modelType(), elements, mapElement);

  addSubElement("V", "transformer_terminal1", Element::STRUCTURE, name(), modelType(), elements, mapElement);
  addSubElement("re", "transformer_terminal1_V", Element::TERMINAL, name(), modelType(), elements, mapElement);
  addSubElement("im", "transformer_terminal1_V", Element::TERMINAL, name(), modelType(), elements, mapElement);

  addSubElement("i", "transformer_terminal1", Element::STRUCTURE, name(), modelType(), elements, mapElement);
  addSubElement("re", "transformer_terminal1_i", Element::TERMINAL, name(), modelType(), elements, mapElement);
  addSubElement("im", "transformer_terminal1_i", Element::TERMINAL, name(), modelType(), elements, mapElement);

  addElement("load", Element::STRUCTURE, elements, mapElement);
  addSubElement("PRefPu", "load", Element::TERMINAL, name(), modelType(), elements, mapElement);
  addSubElement("QRefPu", "load", Element::TERMINAL, name(), modelType(), elements, mapElement);

  addElement("tapChanger", Element::STRUCTURE, elements, mapElement);
  addSubElement("locked", "tapChanger", Element::TERMINAL, name(), modelType(), elements, mapElement);
  addSubElement("deltaUPuTarget", "tapChanger", Element::TERMINAL, name(), modelType(), elements, mapElement);

  addElement("switchOffSignal1", Element::TERMINAL, elements, mapElement);
  addElement("switchOffSignal2", Element::TERMINAL, elements, mapElement);
}

void
ModelLoadOneTransformerTapChanger::defineParameters(vector<ParameterModeler>& parameters) {
  parameters.push_back(ParameterModeler("load_alpha", VAR_TYPE_DOUBLE, EXTERNAL_PARAMETER));
  parameters.push_back(ParameterModeler("load_beta", VAR_TYPE_DOUBLE, EXTERNAL_PARAMETER));
  parameters.push_back(ParameterModeler("transformer_B", VAR_TYPE_DOUBLE, EXTERNAL_PARAMETER));
  parameters.push_back(ParameterModeler("transformer_G", VAR_TYPE_DOUBLE, EXTERNAL_PARAMETER));
  parameters.push_back(ParameterModeler("transformer_X", VAR_TYPE_DOUBLE, EXTERNAL_PARAMETER));
  parameters.push_back(ParameterModeler("transformer_R", VAR_TYPE_DOUBLE, EXTERNAL_PARAMETER));
  parameters.push_back(ParameterModeler("transformer_NbTap", VAR_TYPE_INT, EXTERNAL_PARAMETER));
  parameters.push_back(ParameterModeler("transformer_rTfoMaxPu", VAR_TYPE_DOUBLE, EXTERNAL_PARAMETER));
  parameters.push_back(ParameterModeler("transformer_rTfoMinPu", VAR_TYPE_DOUBLE, EXTERNAL_PARAMETER));
  parameters.push_back(ParameterModeler("transformer_P10Pu", VAR_TYPE_DOUBLE, EXTERNAL_PARAMETER));
  parameters.push_back(ParameterModeler("transformer_Q10Pu", VAR_TYPE_DOUBLE, EXTERNAL_PARAMETER));
  parameters.push_back(ParameterModeler("transformer_U10Pu", VAR_TYPE_DOUBLE, EXTERNAL_PARAMETER));
  parameters.push_back(ParameterModeler("transformer_U1Phase0", VAR_TYPE_DOUBLE, EXTERNAL_PARAMETER));
  parameters.push_back(ParameterModeler("transformer_SNom", VAR_TYPE_DOUBLE, EXTERNAL_PARAMETER));

  parameters.push_back(ParameterModeler("tapChanger_UDeadBand", VAR_TYPE_DOUBLE, EXTERNAL_PARAMETER));
  parameters.push_back(ParameterModeler("tapChanger_UTarget", VAR_TYPE_DOUBLE, EXTERNAL_PARAMETER));
  parameters.push_back(ParameterModeler("tapChanger_t1st", VAR_TYPE_DOUBLE, EXTERNAL_PARAMETER));
  parameters.push_back(ParameterModeler("tapChanger_tNext", VAR_TYPE_DOUBLE, EXTERNAL_PARAMETER));
  parameters.push_back(ParameterModeler("tapChanger_tapMax", VAR_TYPE_INT, EXTERNAL_PARAMETER));
  parameters.push_back(ParameterModeler("tapChanger_tapMin", VAR_TYPE_INT, EXTERNAL_PARAMETER));
  parameters.push_back(ParameterModeler("tapChanger_regulating0", VAR_TYPE_BOOL, EXTERNAL_PARAMETER));
}

void
ModelLoadOneTransformerTapChanger::setSubModelParameters() {
  load_alpha_ = findParameterDynamic("load_alpha").getValue<double>();
  load_half_alpha_ = 0.5 * load_alpha_;
  load_beta_ = findParameterDynamic("load_beta").getValue<double>();
  load_half_beta_ = 0.5 * load_beta_;
  transformer_B_ = findParameterDynamic("transformer_B").getValue<double>();
  transformer_G_ = findParameterDynamic("transformer_G").getValue<double>();
  transformer_X_ = findParameterDynamic("transformer_X").getValue<double>();
  transformer_R_ = findParameterDynamic("transformer_R").getValue<double>();
  transformer_NbTap_ = findParameterDynamic("transformer_NbTap").getValue<int>();
  transformer_rTfoMaxPu_ = findParameterDynamic("transformer_rTfoMaxPu").getValue<double>();
  transformer_rTfoMinPu_ = findParameterDynamic("transformer_rTfoMinPu").getValue<double>();
  transformer_P10Pu_ = findParameterDynamic("transformer_P10Pu").getValue<double>();
  transformer_Q10Pu_ = findParameterDynamic("transformer_Q10Pu").getValue<double>();
  transformer_U10Pu_ = findParameterDynamic("transformer_U10Pu").getValue<double>();
  transformer_U1Phase0_ = findParameterDynamic("transformer_U1Phase0").getValue<double>();
  transformer_SNom_ = findParameterDynamic("transformer_SNom").getValue<double>();

  tapChanger_UDeadBand_ = findParameterDynamic("tapChanger_UDeadBand").getValue<double>();
  tapChanger_UTarget_ = findParameterDynamic("tapChanger_UTarget").getValue<double>();
  tapChanger_t1st_ = findParameterDynamic("tapChanger_t1st").getValue<double>();
  tapChanger_tNext_ = findParameterDynamic("tapChanger_tNext").getValue<double>();
  tapChanger_tapMax_ = findParameterDynamic("tapChanger_tapMax").getValue<int>();
  tapChanger_tapMin_ = findParameterDynamic("tapChanger_tapMin").getValue<int>();
  tapChanger_regulating0_ = findParameterDynamic("tapChanger_regulating0").getValue<bool>();

  tapChanger_increaseTapToIncreaseValue_ = true;
}

void
ModelLoadOneTransformerTapChanger::setFequations() {
  fEquationIndex_[0] = std::string("PPu = Re(V * Conj(i))");
  fEquationIndex_[1] = std::string("QPu = Im(V * Conj(i))");
  fEquationIndex_[2] = std::string("if running PPu = PRefPu * (UPu / U0Pu) ^ alpha) else i_re = 0");
  fEquationIndex_[3] = std::string("if running QPu = QRefPu * (UPu / U0Pu) ^ beta) else i_im = 0");
  fEquationIndex_[4] = std::string("if running Re(i1 = rTfoPu * (YPu * V2 - i2)) else Re(i1 + i2 = 0)");
  fEquationIndex_[5] = std::string("if running Im(i1 = rTfoPu * (YPu * V2 - i2)) else Im(i1 + i2 = 0)");
  fEquationIndex_[6] = std::string("if running Re(ZPu * i1 = rTfoPu * rTfoPu * V1 - rTfoPu * V2) else Re(V1 = V2)");
  fEquationIndex_[7] = std::string("if running Im(ZPu * i1 = rTfoPu * rTfoPu * V1 - rTfoPu * V2) else Im(V1 = V2)");
}

void
ModelLoadOneTransformerTapChanger::setGequations() {
  gEquationIndex_[0] = "((switchoff_signal_1 = true or switchOffSignal2 = true) and running) ||"
                       " ((switchoff_signal_1 = false and switchOffSignal2 = false) && not(running))";
  gEquationIndex_[1] = "tap != pre(tap)";
  gEquationIndex_[2] = "UMonitored > valueMax and running";
  gEquationIndex_[3] = "UMonitored < valueMin and running";
  gEquationIndex_[4] = "moveUp and ((t - whenUp >= t1st and tap == tapWhenUp) or (t - whenLastTap >= tNext and tap != tapWhenUp)) and "
                       "tap < tapMax and not(locked) and running";
  gEquationIndex_[5] = "moveDown and ((t - whenDown >= t1st and tap == tapWhenDown) or (t - whenLastTap >= tNext and tap != tapWhenDown)) and "
                       "tap > tapMin not(locked) and running";
}

}  // namespace DYN
