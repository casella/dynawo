within Dynawo.Examples.ENTSOE;

/*
* Copyright (c) 2015-2021, RTE (http://www.rte-france.com)
* See AUTHORS.txt
* All rights reserved.
* This Source Code Form is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
* file, you can obtain one at http://mozilla.org/MPL/2.0/.
* SPDX-License-Identifier: MPL-2.0
*
* This file is part of Dynawo, an hybrid C++/Modelica open source suite of simulation tools for power systems.
*/

model TestCase3 "Synchronous machine connected to an infinite bus, with governor, PSS and AVR"

  import Modelica;
  import Dynawo;

  extends Modelica.Icons.Example;

  // Grid with impedance
  Dynawo.Electrical.Buses.InfiniteBus infiniteBus(UPhase = 0.000242, UPu = 0.952859) annotation(
    Placement(visible = true, transformation(origin = {-132, 0}, extent = {{-16, -16}, {16, 16}}, rotation = -90)));
  Dynawo.Electrical.Lines.Line gridImpedance(BPu = 0, GPu = 0, RPu = 0.0036, XPu = 0.036) annotation(
    Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  // Transformer
  Dynawo.Electrical.Transformers.TransformerFixedRatio transformer(BPu = 0, GPu = 0, RPu = 0.0003, XPu = 0.032, rTfoPu = 1) annotation(
    Placement(visible = true, transformation(origin = {-32, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  // Generator
  BaseClasses.GeneratorSynchronousInterfaces generatorSynchronous(
   Ce0Pu = 0.95,
   Cm0Pu = 1,
   Cos2Eta0 = 0.586492,
   DPu = 0,
   Efd0Pu = 2.50789,
   ExcitationPu = Dynawo.Electrical.Machines.OmegaRef.BaseClasses.GeneratorSynchronousParameters.ExcitationPuType.NominalStatorVoltageNoLoad,
   H = 4,
   IRotor0Pu = 2.50789,
   IStator0Pu = 5.03993,
   Id0Pu = -0.921348,
   If0Pu = 1.35562,
   Iq0Pu = -0.408844,
   LDPPu = 0.19063,
   LQ1PPu = 0.51659,
   LQ2PPu = 0.24243,
   LambdaAD0Pu = 0.803398,
   LambdaAQ0Pu = -0.674593,
   LambdaAirGap0Pu = 1.04906,
   LambdaD0Pu = 0.803398,
   LambdaQ10Pu = -0.674593,
   LambdaQ20Pu = -0.674593,
   Lambdad0Pu = 0.665196,
   Lambdaf0Pu = 1.10733,
   Lambdaq0Pu = -0.73592,
   LdPPu = 0.15,
   LfPPu = 0.2242,
   LqPPu = 0.15,
   MdPPu = 1.85,
   MdSat0PPu = 1.85,
   Mds0Pu = 1.85,
   Mi0Pu = 1.7673,
   MqPPu = 1.65,
   MqSat0PPu = 1.65,
   Mqs0Pu = 1.65,
   MrcPPu = 0,
   MsalPu = 0.2,
   P0Pu = -4.75,
   PGen0Pu = 4.75,
   PNomAlt = 475,
   PNomTurb = 475,
   Pm0Pu = 1,
   Q0Pu = -1.56,
   QGen0Pu = 1.56,
   QStator0Pu = 1.56,
   RDPPu = 0.02933,
   RQ1PPu = 0.0035,
   RQ2PPu = 0.02227,
   RTfPu = 0,
   RaPPu = 0,
   RfPPu = 0.00128,
   SNom = 500,
   Sin2Eta0 = 0.413508,
   SnTfo = 500,
   Theta0 = 0.996978,
   ThetaInternal0 = 0.996978,
   U0Pu = 0.992,
   UBaseHV = 400,
   UBaseLV = 21,
   UNom = 21,
   UNomHV = 400,
   UNomLV = 21,
   UPhase0 = 0.161146,
   UStator0Pu = 0.992,
   Ud0Pu = 0.73592,
   Uf0Pu = 0.00173519,
   Uq0Pu = 0.665196,
   XTfPu = 0,
   md = 0,
   mq = 0,
   nd = 0,
   nq = 0) annotation(
    Placement(visible = true, transformation(origin = {20, 1.9984e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Dynawo.Electrical.Controls.Basics.SetPoint Omega0Pu(Value0 = 1);

  // Load
  Dynawo.Electrical.Loads.LoadAlphaBeta load(alpha = 2, beta = 2, u0Pu = Complex(0.952267, 0)) annotation(
    Placement(visible = true, transformation(origin = {-80, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Dynawo.Electrical.Controls.Basics.SetPoint PRefPu(Value0 = 4.75);
  Dynawo.Electrical.Controls.Basics.SetPoint QRefPu(Value0 = 0.76);
  Dynawo.Electrical.Events.NodeFault nodeFault(RPu = 0.000173, XPu = 0, tBegin = 0.1, tEnd = 0.2) annotation(
    Placement(visible = true, transformation(origin = {-52, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

  // Regulations
  Dynawo.Electrical.Controls.Machines.Governors.IEEE_TGOV1 governor(Dt = 0, Pm0Pu = generatorSynchronous.Pm0Pu, R = 0.05, Tg1 = 0.5, Tg2 = 3, Tg3 = 10, VMax = 1, VMin = 0) annotation(
    Placement(visible = true, transformation(origin = {90, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Dynawo.Electrical.Controls.Machines.PowerSystemStabilizers.IEEE_PSS2A pss(IC1 = 1, IC2 = 3, Ks1 = 10, Ks2 = 0.1564, Ks3 = 1, PGen0Pu = -generatorSynchronous.P0Pu, PNomAlt = generatorSynchronous.PNomAlt, T1 = 0.25, T2 = 0.03, T3 = 0.15, T4 = 0.015, T6 = 1e-5, T7 = 2, T8 = 0.5, T9 = 0.1, Tw1 = 2, Tw2 = 2, Tw3 = 2, Tw4 = 1e-5, Upss0Pu = 0, VstMax = 0.1, VstMin = -0.1) annotation(
    Placement(visible = true, transformation(origin = {90, 0}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Dynawo.Electrical.Controls.Machines.VoltageRegulators.IEEE_SEXS avr(EMax = 4, EMin = 0, Efd0Pu = generatorSynchronous.Efd0Pu, K = 200, Ta = 3, Tb = 10, Te = 0.05, Upss0Pu = 0, Us0Pu = 0.992) annotation(
    Placement(visible = true, transformation(origin = {130, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const(k = avr.UsRef0Pu) annotation(
    Placement(visible = true, transformation(origin = {10, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

equation
  connect(transformer.terminal2, generatorSynchronous.terminal) annotation(
    Line(points = {{-12, 0}, {20, 0}}, color = {0, 0, 255}));
  connect(gridImpedance.terminal2, transformer.terminal1) annotation(
    Line(points = {{-80, 0}, {-52, 0}}, color = {0, 0, 255}));
  connect(gridImpedance.terminal1, infiniteBus.terminal) annotation(
    Line(points = {{-120, 0}, {-132, 0}}, color = {0, 0, 255}));
  connect(load.terminal, gridImpedance.terminal2) annotation(
    Line(points = {{-80, -38}, {-80, 0}}, color = {0, 0, 255}));
  connect(generatorSynchronous.omegaRefPu, Omega0Pu.setPoint) annotation(
    Line);
  connect(load.PRefPu, PRefPu.setPoint) annotation(
    Line);
  connect(load.QRefPu, QRefPu.setPoint) annotation(
    Line);
  connect(nodeFault.terminal, transformer.terminal1) annotation(
    Line(points = {{-52, 50}, {-52, 0}}, color = {0, 0, 255}));
  gridImpedance.switchOffSignal1.value = false;
  gridImpedance.switchOffSignal2.value = false;
  transformer.switchOffSignal1.value = false;
  transformer.switchOffSignal2.value = false;
  load.switchOffSignal1.value = false;
  load.switchOffSignal2.value = false;
  generatorSynchronous.switchOffSignal1.value = false;
  generatorSynchronous.switchOffSignal2.value = false;
  generatorSynchronous.switchOffSignal3.value = false;
  connect(generatorSynchronous.omegaPu_out, governor.omegaPu) annotation(
    Line(points = {{38, 6}, {60, 6}, {60, -36}, {78, -36}}, color = {0, 0, 127}));
  connect(generatorSynchronous.PGenPu_out, pss.PGenPu) annotation(
    Line(points = {{38, 10}, {70, 10}, {70, 6}, {78, 6}}, color = {0, 0, 127}));
  connect(generatorSynchronous.omegaPu_out, pss.omegaPu) annotation(
    Line(points = {{38, -6}, {78, -6}}, color = {0, 0, 127}));
  connect(pss.UpssPu, avr.UpssPu) annotation(
    Line(points = {{101, 0}, {110, 0}, {110, 12}, {118, 12}}, color = {0, 0, 127}));
  connect(generatorSynchronous.UsPu_out, avr.UsPu) annotation(
    Line(points = {{38, 18}, {118, 18}}, color = {0, 0, 127}));
  connect(const.y, avr.UsRefPu) annotation(
    Line(points = {{21, 60}, {80, 60}, {80, 24}, {118, 24}}, color = {0, 0, 127}));
  connect(governor.PmPu, generatorSynchronous.PmPu_in) annotation(
    Line(points = {{102, -30}, {102, -30.5}, {110, -30.5}, {110, -51}, {32, -51}, {32, -16}}, color = {0, 0, 127}));
  connect(avr.EfdPu, generatorSynchronous.efdPu_in) annotation(
    Line(points = {{141, 18}, {150, 18}, {150, -60}, {8, -60}, {8, -16}}, color = {0, 0, 127}));

  annotation(
    experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-06),
    __OpenModelica_simulationFlags(initialStepSize = "0.001", lv = "LOG_STATS", nls = "kinsol", s = "ida", nlsLS = "klu", maxIntegrationOrder = "2", maxStepSize = "10", emit_protected = "()"),
  Diagram(coordinateSystem(extent = {{-160, -100}, {160, 100}})),
  Icon);
end TestCase3;
