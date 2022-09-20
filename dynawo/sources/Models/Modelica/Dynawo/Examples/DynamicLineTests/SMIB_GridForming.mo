within Dynawo.Examples.DynamicLineTests;

/*
* Copyright (c) 2022, RTE (http://www.rte-france.com)
* See AUTHORS.txt
* All rights reserved.
* This Source Code Form is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
* file, you can obtain one at http://mozilla.org/MPL/2.0/.
* SPDX-License-Identifier: MPL-2.0
*
* This file is part of Dynawo, an hybrid C++/Modelica open source suite
* of simulation tools for power systems.
*/

model SMIB_GridForming "Node fault on a line for equivalent SMIB model using a synchronous machine and a grid forming inverter"

  import Dynawo;
  import Modelica;

  extends Icons.Example;

  parameter Real x = 0.5 "Emplacement of the fault relative to the line lenght : x= default location /line lenght";
  parameter Types.PerUnit UBusPu = 0.604201 "Complex voltage amplitude at the infinite bus in pu (base Unom)";
  parameter Types.Angle UBusPhase = 0.022836 "Complex voltage phase at the infinite bus in rad";
  parameter Types.PerUnit RPu = 0.00375 "Resistance in pu (base SnRef, UNom) ";
  parameter Types.PerUnit XPu = 0.0375 "Reactance in pu (base SnRef, UNom)";
  parameter Types.PerUnit GPu = 0 "Half-conductance in pu (base SnRef, UNom)";
  parameter Types.PerUnit BPu = 0 "Half-susceptance in pu (base SnRef, UNom)";
  parameter Real delta_t = 0.4 "Fault duration in seconds";
  parameter Real coeff_machine = 0.4 "Synchronous machine percentage (between 0 and 1)";
  parameter Real conv = 0.6 "Converter percentage (between 0 and 1)";
  parameter Real teta = 0.600788 "Start value of the phase shift between the converter's rotating frame and the grid rotating frame (in rad)";
  parameter Real P_mec = 0.9030004 "Start value of mechanical power in pu (base PNomTurb/OmegaNom)";
  parameter Real efd = 2.465967 "Start value of input exciter voltage in pu (user-selected base)";
  parameter Real Udc0 = 1.088987 "DC Voltage reference in pu (base UNom)";
  parameter Real Idc0 = 0.83611 "DC Current reference in pu (base UNom)";
  parameter Real Uref0 = 1.088987 "Reference DC voltage on the DC side in pu (base UNom, SNom)";
  parameter Real PRef250 = 0.910513 "Active power reference at the converter's capacitor in pu (base SNom)";
  parameter Real QRef250 = 0.6462967 "Reactive power reference at the converter's capacitor in pu (base SNom)";
  parameter Real omegaRefPu = 1 "Grid frequency in pu";


  //Generator, lines and infinite bus
  Dynawo.Examples.BaseClasses.GeneratorSynchronousInterfaces generatorSynchronous(Ce0Pu = 0.9030004, Cm0Pu = 0.9030004, Cos2Eta0 = 0.6889078, DPu = 0, Efd0Pu = efd, ExcitationPu = Dynawo.Electrical.Machines.OmegaRef.BaseClasses.GeneratorSynchronousParameters.ExcitationPuType.NominalStatorVoltageNoLoad, H = 3.5, IRotor0Pu = 2.465967, IStator0Pu = 8.880566, Id0Pu = -0.9197753, If0Pu = 1.485522, Iq0Pu = -0.3926077, LDPPu = 0.16634, LQ1PPu = 0.92815, LQ2PPu = 0.12046, LambdaAD0Pu = 0.8934501, LambdaAQ0Pu = -0.6003912, LambdaAirGap0Pu = 1.07644, LambdaD0Pu = 0.8934501, LambdaQ10Pu = -0.6003912, LambdaQ20Pu = -0.6003912, Lambdad0Pu = 0.7554838, Lambdaf0Pu = 1.14584, Lambdaq0Pu = -0.6592823, LdPPu = 0.15, LfPPu = 0.16990, LqPPu = 0.15, MdPPu = 1.66, MdSat0PPu = 1.579239, Mds0Pu = 1.578475, Mi0Pu = 1.563685, MqPPu = 1.61, MqSat0PPu = 1.529239, Mqs0Pu = 1.530931, MrcPPu = 0, MsalPu = 0.05, P0Pu = -19.98 * coeff_machine , PGen0Pu = 19.98 * coeff_machine, PNomAlt = 2200 * coeff_machine, PNomTurb = 2220*coeff_machine, Pm0Pu = P_mec, Q0Pu = -9.68 * coeff_machine, QGen0Pu = 9.6789 * coeff_machine, QStator0Pu = 3.872, RDPPu = 0.03339, RQ1PPu = 0.00924, RQ2PPu = 0.02821, RTfPu = 0, RaPPu = 0.003, RfPPu = 0.00074, SNom = 2220*coeff_machine , Sin2Eta0 = 0.3110922, SnTfo = 2220*coeff_machine, Theta0 = 1.2062, ThetaInternal0 = 0.7161999 , U0Pu = 1, UBaseHV = 24, UBaseLV = 24, UNom = 24, UNomHV = 24, UNomLV = 24, UPhase0 = 0.494442, UStator0Pu = 1.0, Ud0Pu = 0.656523, Uf0Pu = 0.001099287, Uq0Pu = 0.754306, XTfPu = 0, md = 0.031, mq = 0.031, nd = 6.93, nq = 6.93) annotation(
    Placement(visible = true, transformation(origin = {74, 66}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Dynawo.Electrical.Transformers.TransformerFixedRatio transformer(BPu = 0, GPu = 0, RPu = 0, XPu = 0.00675, rTfoPu = 1) annotation(
    Placement(visible = true, transformation(origin = {30, 66}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
  Dynawo.Electrical.Lines.Line line2(BPu = BPu, GPu = GPu, RPu = RPu, XPu = XPu) annotation(
    Placement(visible = true, transformation(origin = {-32, 44}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
  Dynawo.Electrical.Lines.Line line1(BPu = BPu * (1 - x), GPu = GPu, RPu = RPu * (1 - x), XPu = XPu * (1 - x)) annotation(
    Placement(visible = true, transformation(origin = {-68, 90}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
  Dynawo.Electrical.Controls.Basics.Step PmPu(Value0 = P_mec, Height = 0.02, tStep = 1000);
  Dynawo.Electrical.Controls.Basics.SetPoint Omega0Pu(Value0 = 1);
  Dynawo.Electrical.Controls.Basics.SetPoint EfdPu(Value0 = efd);
  Dynawo.Electrical.Lines.Line line(BPu = BPu * x, GPu = GPu, RPu = RPu * x, XPu = XPu * x) annotation(
    Placement(visible = true, transformation(origin = {-8.88178e-16, 90}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
  Dynawo.Electrical.Buses.InfiniteBus infiniteBus(UPhase = UBusPhase, UPu = UBusPu) annotation(
    Placement(visible = true, transformation(origin = {-84, 68}, extent = {{-8, -8}, {8, 8}}, rotation = -90)));
  Dynawo.Electrical.Events.NodeFault nodeFault(RPu = 0.00001, XPu = 0.00001, tBegin = 5, tEnd = 5 + delta_t) annotation(
    Placement(visible = true, transformation(origin = {-42, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Step UFilterRef250Pu(height = 0, offset = Uref0, startTime = 1) annotation(
    Placement(visible = true, transformation(origin = {-77, -63}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  Modelica.Blocks.Sources.Step IdcSourceRef250Pu(height = 0, offset = Idc0, startTime = 1) annotation(
    Placement(visible = true, transformation(origin = {-67, -82}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  Dynawo.Electrical.Controls.Converters.GridFormingControlDroopControl Droop(Cfilter = 0.066, IdConv0Pu = 0.83611, IdPcc0Pu = 0.83611, IdcSource0Pu = 0.840569, IdcSourceRef0Pu = Idc0, Imax = 1, IqConv0Pu = -0.5216111, IqPcc0Pu = -0.5934843, Kff = 0.01, Kic = 1.19, Kiv = 1.161022, KpVI = 0.67, Kpc = 0.7388, Kpdc = 50, Kpv = 0.52, Lfilter = 0.15, Mp = 0.01 , Mq = 0.01, PRef0Pu = PRef250, QRef0Pu = QRef250, Rfilter = 0.005, Theta0 = teta, UdConv0Pu = 1.171409, UdFilter0Pu = Uref0, UdcSource0Pu = Udc0, UdcSourcePu(fixed = false), UqConv0Pu = 0.1228084, Wf = 60, Wff = 16.66, XRratio = 5, currentLoop(integratord(y_start = 0.00323126), integratorq(y_start = -0.000164394)), droopControl(firstOrder(y_start = -7.3445e-5), firstOrder1(y_start = 0.102988), firstOrder2(y_start = 0.00622874), firstOrder3(y_start = -0.0010158), integrator(y_start = -0.0502873)), idConvPu(fixed = false), idPccPu(fixed = false), iqConvPu(fixed = false), iqPccPu(fixed = false), udFilterPu(fixed = false), uqFilterPu(fixed = false)) annotation(
    Placement(visible = true, transformation(origin = {-18, -68}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
  Modelica.Blocks.Sources.Step UdcSourceRef250Pu(height = 0, offset = Udc0, startTime = 1) annotation(
    Placement(visible = true, transformation(origin = {-49, -95}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  Modelica.Blocks.Sources.Step PRef250Pu(height = 0, offset = PRef250, startTime = 1) annotation(
    Placement(visible = true, transformation(origin = {-77, -23}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  Dynawo.Electrical.Sources.Converter Conv250(Cdc = 0.01, Cfilter = 0.066, IdConv0Pu = 0.83611, IdPcc0Pu = 0.83611, IdcSource0Pu = 0.840569, IqConv0Pu = -0.5216111, IqPcc0Pu = -0.5934843, Lfilter = 0.15, Ltransformer = 0.2, PGenPu(fixed = true, start = 19.98*conv), QGenPu(fixed = true, start = 9.68 * conv ), Rfilter = 0.005, Rtransformer = 0.01, SNom = 2220*conv, Theta0 = teta, UdConv0Pu = 1.171409, UdFilter0Pu = 1.088987, UdPcc0Pu = 0.9619292, UdcSource0Pu = 1.088987, UqConv0Pu = 0.1228084, UqPcc0Pu = -0.1612872, i0Pu = Complex(-13.65555, 0.2252762), theta(fixed = true, start = teta ), u0Pu = Complex(0.8846606, 0.4107273)) annotation(
    Placement(visible = true, transformation(origin = {28, -68}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
  Modelica.Blocks.Sources.Step QRef250Pu(height = 0, offset = QRef250, startTime = 1) annotation(
    Placement(visible = true, transformation(origin = {-77, -42}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  Dynawo.Electrical.Controls.Machines.Governors.Simplified.GoverProportional goverProportional(KGover = 1, PMax = 3000, PMin = 0, PNom = 1000, Pm0Pu = generatorSynchronous.Pm0Pu) annotation(
    Placement(visible = true, transformation(origin = {43, 19}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  Dynawo.Electrical.Controls.Machines.VoltageRegulators.Simplified.VRProportionalIntegral vRProportionalIntegral(Efd0Pu = generatorSynchronous.Efd0Pu, EfdMaxPu = 5, EfdMinPu = -5, Gain = 20, LagEfdMax = 0, LagEfdMin = 0, Us0Pu = generatorSynchronous.U0Pu, UsRef0Pu = generatorSynchronous.U0Pu, UsRefMaxPu = 5, UsRefMinPu = -5, tIntegral = 1, yIntegrator0 = generatorSynchronous.Efd0Pu) annotation(
    Placement(visible = true, transformation(origin = {43, 3}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant PmRefPu(k = generatorSynchronous.Pm0Pu) annotation(
    Placement(visible = true, transformation(origin = {-24, 28}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant UsRefPu(k = generatorSynchronous.U0Pu) annotation(
    Placement(visible = true, transformation(origin = {-24, -4}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));

initial equation
  der(generatorSynchronous.lambdafPu) = 0;
  der(generatorSynchronous.lambdaDPu) = 0;
  der(generatorSynchronous.lambdaQ1Pu) = 0;
  der(generatorSynchronous.lambdaQ2Pu) = 0;
  der(generatorSynchronous.theta) = 0;
  der(Conv250.theta) = 0;
  der(generatorSynchronous.omegaPu.value) = 0;

equation
  connect(transformer.terminal2, generatorSynchronous.terminal) annotation(
    Line(points = {{44, 66}, {74, 66}}, color = {0, 0, 255}));
  generatorSynchronous.omegaRefPu.value = omegaRefPu;
  Droop.omegaRefPu = omegaRefPu;
  line1.switchOffSignal1.value = false;
  Conv250.switchOffSignal1.value = false;
  Conv250.switchOffSignal2.value = false;
  Conv250.switchOffSignal3.value = false;
  line1.switchOffSignal2.value = false;
  line2.switchOffSignal1.value = false;
  line2.switchOffSignal2.value = false;
  line.switchOffSignal1.value = false;
  line.switchOffSignal2.value = false;
  transformer.switchOffSignal1.value = false;
  transformer.switchOffSignal2.value = false;
  generatorSynchronous.switchOffSignal1.value = false;
  generatorSynchronous.switchOffSignal2.value = false;
  generatorSynchronous.switchOffSignal3.value = false;
  connect(line2.terminal2, transformer.terminal1) annotation(
    Line(points = {{-16, 44}, {16, 44}, {16, 66}}, color = {0, 0, 255}));
  connect(line1.terminal2, line.terminal1) annotation(
    Line(points = {{-52, 90}, {-16, 90}}, color = {0, 0, 255}));
  connect(line.terminal2, transformer.terminal1) annotation(
    Line(points = {{16, 90}, {16, 66}}, color = {0, 0, 255}));
  connect(line1.terminal1, infiniteBus.terminal) annotation(
    Line(points = {{-84, 90}, {-84, 68}}, color = {0, 0, 255}));
  connect(line2.terminal1, infiniteBus.terminal) annotation(
    Line(points = {{-48, 44}, {-84, 44}, {-84, 68}}, color = {0, 0, 255}));
  connect(nodeFault.terminal, line.terminal1) annotation(
    Line(points = {{-42, 66}, {-42, 83}, {-16, 83}, {-16, 90}}, color = {0, 0, 255}));
  connect(Droop.UdcSourceRefOutPu, Conv250.UdcSourceRefPu) annotation(
    Line(points = {{-3.3, -77.8}, {12.7, -77.8}}, color = {0, 0, 127}));
  connect(Droop.uqConvRefPu, Conv250.uqConvRefPu) annotation(
    Line(points = {{-3.3, -73.6}, {12.7, -73.6}}, color = {0, 0, 127}));
  connect(Conv250.UdcSourcePu, Droop.UdcSourcePu) annotation(
    Line(points = {{42.7, -68}, {68.85, -68}, {68.85, -37}, {-18.3, -37}, {-18.3, -53}}, color = {0, 0, 127}));
  connect(UFilterRef250Pu.y, Droop.UFilterRefPu) annotation(
    Line(points = {{-71.5, -63}, {-55.25, -63}, {-55.25, -68}, {-33, -68}}, color = {0, 0, 127}));
  connect(UdcSourceRef250Pu.y, Droop.UdcSourceRefPu) annotation(
    Line(points = {{-43.5, -95}, {-37.5, -95}, {-37.5, -82}, {-33, -82}}, color = {0, 0, 127}));
  connect(IdcSourceRef250Pu.y, Droop.IdcSourceRefPu) annotation(
    Line(points = {{-61.5, -82}, {-50, -82}, {-50, -76}, {-33, -76}}, color = {0, 0, 127}));
  connect(Droop.idPccPu, Conv250.idPccPu) annotation(
    Line(points = {{-9.6, -53.3}, {-9.6, -46.6}, {60.8, -46.6}, {60.8, -60.3}, {43.4, -60.3}}, color = {0, 0, 127}));
  connect(Conv250.iqPccPu, Droop.iqPccPu) annotation(
    Line(points = {{42.7, -76.4}, {78.85, -76.4}, {78.85, -28.2}, {-26.3, -28.2}, {-26.3, -53.4}}, color = {0, 0, 127}));
  connect(PRef250Pu.y, Droop.PRefPu) annotation(
    Line(points = {{-71.5, -23}, {-44, -23}, {-44, -54}, {-33, -54}}, color = {0, 0, 127}));
  connect(Droop.omegaPu, Conv250.omegaPu) annotation(
    Line(points = {{-3.3, -80.6}, {12.7, -80.6}}, color = {0, 0, 127}));
  connect(QRef250Pu.y, Droop.QRefPu) annotation(
    Line(points = {{-71.5, -42}, {-50, -42}, {-50, -60}, {-33, -60}}, color = {0, 0, 127}));
  connect(Droop.theta, Conv250.theta) annotation(
    Line(points = {{-3.3, -55.4}, {12.7, -55.4}}, color = {0, 0, 127}));
  connect(Droop.idConvPu, Conv250.idConvPu) annotation(
    Line(points = {{-13.8, -53.3}, {-13.8, -42.6}, {64.4, -42.6}, {64.4, -64.3}, {43.2, -64.3}}, color = {0, 0, 127}));
  connect(Droop.udFilterPu, Conv250.udFilterPu) annotation(
    Line(points = {{-5.82, -53.3}, {-5.82, -50.6}, {56.36, -50.6}, {56.36, -55.3}, {43.18, -55.3}}, color = {0, 0, 127}));
  connect(Droop.IdcSourcePu, Conv250.IdcSourcePu) annotation(
    Line(points = {{-3.3, -68}, {12.7, -68}}, color = {0, 0, 127}));
  connect(Conv250.uqFilterPu, Droop.uqFilterPu) annotation(
    Line(points = {{42.7, -80.6}, {83.85, -80.6}, {83.85, -20.8}, {-31.3, -20.8}, {-31.3, -52.6}}, color = {0, 0, 127}));
  connect(Conv250.iqConvPu, Droop.iqConvPu) annotation(
    Line(points = {{42.7, -72.2}, {73.85, -72.2}, {73.85, -32.6}, {-23.3, -32.6}, {-23.3, -53.2}}, color = {0, 0, 127}));
  connect(Droop.udConvRefPu, Conv250.udConvRefPu) annotation(
    Line(points = {{-3.3, -62.4}, {12.7, -62.4}}, color = {0, 0, 127}));
  connect(Conv250.terminal, transformer.terminal1) annotation(
    Line(points = {{28, -83}, {16, -83}, {16, 66}}, color = {0, 0, 255}));
  connect(UsRefPu.y, vRProportionalIntegral.UsRefPu) annotation(
    Line(points = {{-17, -4}, {10.5, -4}, {10.5, 3}, {36, 3}}, color = {0, 0, 127}));
  connect(goverProportional.PmRefPu, PmRefPu.y) annotation(
    Line(points = {{35.5, 22}, {9.25, 22}, {9.25, 28}, {-17, 28}}, color = {0, 0, 127}));
  connect(goverProportional.PmPu, generatorSynchronous.PmPu_in) annotation(
    Line(points = {{49, 19}, {86, 19}, {86, 50}}, color = {0, 0, 127}));
  connect(vRProportionalIntegral.EfdPu, generatorSynchronous.efdPu_in) annotation(
    Line(points = {{50, 4}, {62, 4}, {62, 50}}, color = {0, 0, 127}));
  connect(generatorSynchronous.omegaPu_out, goverProportional.omegaPu) annotation(
    Line(points = {{92, 60}, {98, 60}, {98, 10}, {30, 10}, {30, 16}, {35, 16}}, color = {0, 0, 127}));
  connect(generatorSynchronous.UsPu_out, vRProportionalIntegral.UsPu) annotation(
    Line(points = {{92, 84}, {96, 84}, {96, -4}, {26, -4}, {26, 0}, {40, 0}}, color = {0, 0, 127}));


  annotation(
    preferredView = "text",
    Documentation(info = "<html><head></head><body>
The purpose of this test case is to evaluate the transient stability of an equivalent SMIB model using a synchronous machine and a grid forming inverter by evaluating the rotor angle Theta for a node fault starting at (t_begin = 1 s) and (t_end = t_begin + delta_t). The following figures show the excepted evolution of the generator rotor angle a if delta_t < CCT (the critical clearing time).

    <figure>
    <img width=\"350\" src=\"modelica://Dynawo/Examples/DynamicLineTests/Images/theta_phasor.png\">
    </figure>



For delta_t > CCT, we represent the loss of synchronism with an assert stopping the simulation if theta exceeds 270 degrees.
</body></html>"),
    experiment(StartTime = 0, StopTime = 20, Tolerance = 1e-06, Interval = 0.04));
end SMIB_GridForming;
