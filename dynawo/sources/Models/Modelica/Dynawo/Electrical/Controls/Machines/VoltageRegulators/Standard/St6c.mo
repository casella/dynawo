within Dynawo.Electrical.Controls.Machines.VoltageRegulators.Standard;

/*
* Copyright (c) 2024, RTE (http://www.rte-france.com)
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

model St6c "IEEE exciter type ST6C model"

  //Regulation parameters
  parameter Types.CurrentModulePu IlrPu "Rotor current threshold of field current limiter in pu (base SNom, user-selected base voltage)";
  parameter Types.PerUnit Kc "Rectifier loading factor proportional to commutating reactance";
  parameter Types.PerUnit Kcl "Field current limiter conversion factor";
  parameter Types.PerUnit Kff "Feedforward gain of inner loop field regulator";
  parameter Types.PerUnit Kg "Feedback gain constant of inner loop field regulator";
  parameter Types.PerUnit Ki "Potential circuit (current) gain coefficient";
  parameter Types.PerUnit Kia "Integral gain of PI";
  parameter Types.PerUnit Klr "Gain of field current limiter";
  parameter Types.PerUnit Km "Gain of error of inner loop field regulator";
  parameter Types.PerUnit Kp "Potential circuit gain";
  parameter Types.PerUnit Kpa "Proportional gain of PI";
  parameter Integer PositionOel "Input location : (0) none, (1) voltage error summation, (2) take-over at AVR input, (3) AVR input summation, (4) take-over at AVR output";
  parameter Integer PositionScl "Input location : (0) none, (1) voltage error summation, (2) take-over at AVR input, (3) AVR input summation, (4) take-over at AVR output";
  parameter Integer PositionUel "Input location : (0) none, (1) voltage error summation, (2) take-over at AVR input, (3) AVR input summation, (4) take-over at AVR output";
  parameter Boolean Sw1 "If true, power source derived from terminal voltage, if false, independent from terminal voltage";
  parameter Types.Time tA "Voltage regulator time constant in s";
  parameter Types.Time tG "Feedback time constant of inner loop field regulator in s";
  parameter Types.Angle Thetap "Potential circuit phase angle in rad";
  parameter Types.Time tR "Stator voltage filter time constant in s";
  parameter Types.VoltageModulePu VaMaxPu "Maximum output voltage of limited first order in pu";
  parameter Types.VoltageModulePu VaMinPu "Minimum output voltage of limited first order in pu";
  parameter Types.VoltageModulePu VbMaxPu "Maximum available exciter field voltage in pu (base UNom)";
  parameter Types.VoltageModulePu VmMaxPu "Maximum output voltage of second PI in pu";
  parameter Types.VoltageModulePu VmMinPu "Minimum output voltage of second PI in pu";
  parameter Types.VoltageModulePu VrMaxPu "Maximum output voltage of first PI in pu";
  parameter Types.VoltageModulePu VrMinPu "Minimum output voltage of first PI in pu";
  parameter Types.PerUnit XlPu "Reactance associated with potential source in pu (base SNom, UNom)";

  //Input variables
  Modelica.Blocks.Interfaces.RealInput IrPu(start = Ir0Pu) "Rotor current in pu (base SNom, user-selected base voltage)" annotation(
    Placement(visible = true, transformation(origin = {-460, 200}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, 100}, extent = {{20, -20}, {-20, 20}}, rotation = 180)));
  Modelica.ComplexBlocks.Interfaces.ComplexInput itPu(re(start = it0Pu.re), im(start = it0Pu.im)) "Complex stator current in pu (base UNom)" annotation(
    Placement(visible = true, transformation(origin = {-460, 120}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {120, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 180)));
  Modelica.Blocks.Interfaces.RealInput UOelPu(start = UOel0Pu) "Overexcitation limitation output voltage in pu (base UNom)" annotation(
    Placement(visible = true, transformation(origin = {-460, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput UPssPu(start = 0) "Power system stabilizer output voltage in pu (base UNom)" annotation(
    Placement(visible = true, transformation(origin = {-460, -200}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, -100}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput USclOelPu(start = USclOel0Pu) "Stator current overexcitation limitation output voltage in pu (base UNom)" annotation(
    Placement(visible = true, transformation(origin = {-460, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {120, 80}, extent = {{20, -20}, {-20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput USclUelPu(start = USclUel0Pu) "Stator current underexcitation limitation output voltage in pu (base UNom)" annotation(
    Placement(visible = true, transformation(origin = {-460, -160}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {120, 40}, extent = {{20, -20}, {-20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput UsPu(start = Us0Pu) "Stator voltage in pu (base UNom)" annotation(
    Placement(visible = true, transformation(origin = {-460, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, -20}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput UsRefPu(start = Us0Pu) "Control voltage in pu (base UNom)" annotation(
    Placement(visible = true, transformation(origin = {-460, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, 20}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.ComplexBlocks.Interfaces.ComplexInput utPu(re(start = ut0Pu.re), im(start = ut0Pu.im)) "Complex stator voltage in pu (base UNom)" annotation(
    Placement(visible = true, transformation(origin = {-460, 160}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {120, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 180)));
  Modelica.Blocks.Interfaces.RealInput UUelPu(start = UUel0Pu) "Underexcitation limitation output voltage in pu (base UNom)" annotation(
    Placement(visible = true, transformation(origin = {-460, -120}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  //Output variable
  Modelica.Blocks.Interfaces.RealOutput EfdPu(start = Efd0Pu) "Excitation voltage in pu (user-selected base voltage)" annotation(
    Placement(visible = true, transformation(origin = {430, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

  Modelica.Blocks.Sources.Constant const(k = VbMaxPu) annotation(
    Placement(visible = true, transformation(origin = {-70, 200}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Division division annotation(
    Placement(visible = true, transformation(origin = {-170, 160}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Min min5 annotation(
    Placement(visible = true, transformation(origin = {-10, 160}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Product product1 annotation(
    Placement(visible = true, transformation(origin = {-70, 140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Dynawo.Electrical.Controls.Machines.VoltageRegulators.Standard.BaseClasses.RectifierRegulationCharacteristic rectifierRegulationCharacteristic annotation(
    Placement(visible = true, transformation(origin = {-130, 160}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Dynawo.NonElectrical.Blocks.Continuous.PIAntiWindupVariableLimits limPI1(Ki = Kia, Kp = Kpa, Y0 = Kg * Efd0Pu) annotation(
    Placement(visible = true, transformation(origin = {-50, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.FirstOrder firstOrder(T = tR, y_start = Us0Pu) annotation(
    Placement(visible = true, transformation(origin = {-410, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Product product annotation(
    Placement(visible = true, transformation(origin = {370, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Feedback feedback annotation(
    Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain1(k = Kc) annotation(
    Placement(visible = true, transformation(origin = {-230, 200}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Dynawo.NonElectrical.Blocks.NonLinear.LimitedFirstOrder limitedFirstOrder(Y0 = Efd0Pu / Vb0Pu, YMax = VmMaxPu, YMin = VmMinPu, tFilter = tA) annotation(
    Placement(visible = true, transformation(origin = {310, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add add(k2 = -1) annotation(
    Placement(visible = true, transformation(origin = {-350, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add add2(k1 = Kff, k2 = Km) annotation(
    Placement(visible = true, transformation(origin = {50, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.Limiter limiter(uMax = VrMaxPu, uMin = VrMinPu) annotation(
    Placement(visible = true, transformation(origin = {90, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Min min4 annotation(
    Placement(visible = true, transformation(origin = {270, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.FirstOrder firstOrder1(T = tG, k = Kg, y_start = Efd0Pu) annotation(
    Placement(visible = true, transformation(origin = {310, -200}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Dynawo.Electrical.Controls.Machines.VoltageRegulators.Standard.BaseClasses.PotentialCircuit potentialCircuit(Ki = Ki, Kp = Kp, Theta = Thetap, X = XlPu) annotation(
    Placement(visible = true, transformation(origin = {-390, 140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Feedback feedback1 annotation(
    Placement(visible = true, transformation(origin = {-280, 20}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const1(k = IlrPu * Kcl) annotation(
    Placement(visible = true, transformation(origin = {-330, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain(k = Klr) annotation(
    Placement(visible = true, transformation(origin = {-230, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Max max3 annotation(
    Placement(visible = true, transformation(origin = {-170, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const2(k = VrMinPu) annotation(
    Placement(visible = true, transformation(origin = {-230, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Min min3 annotation(
    Placement(visible = true, transformation(origin = {-110, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const3(k = VaMaxPu) annotation(
    Placement(visible = true, transformation(origin = {-170, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const4(k = VaMinPu) annotation(
    Placement(visible = true, transformation(origin = {-110, -200}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Sum sum1(nin = 4) annotation(
    Placement(visible = true, transformation(origin = {-290, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Sum sum2(nin = 5) annotation(
    Placement(visible = true, transformation(origin = {-110, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.MinMax max1(nu = 3) annotation(
    Placement(visible = true, transformation(origin = {-230, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.MinMax min1(nu = 3) annotation(
    Placement(visible = true, transformation(origin = {-170, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.MinMax max2(nu = 3) annotation(
    Placement(visible = true, transformation(origin = {150, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.MinMax min2(nu = 3) annotation(
    Placement(visible = true, transformation(origin = {210, -74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.Switch switch annotation(
    Placement(visible = true, transformation(origin = {-330, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const5(k = Kp) annotation(
    Placement(visible = true, transformation(origin = {-390, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanConstant booleanConstant(k = Sw1) annotation(
    Placement(visible = true, transformation(origin = {-390, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

  //Generator initial parameters
  parameter Types.VoltageModulePu Efd0Pu "Initial excitation voltage in pu (user-selected base voltage)";
  parameter Types.CurrentModulePu Ir0Pu "Initial rotor current in pu (base SNom, user-selected base voltage)";
  parameter Types.ComplexCurrentPu it0Pu "Initial complex stator current in pu (base SNom, UNom)";
  parameter Types.VoltageModulePu Us0Pu "Initial stator voltage in pu (base UNom)";
  parameter Types.ComplexVoltagePu ut0Pu "Initial complex stator voltage in pu (base UNom)";

  //Initial parameters (inputs)
  parameter Types.VoltageModulePu UOel0Pu = 0 "Overexcitation limitation initial output voltage in pu (base UNom)";
  parameter Types.VoltageModulePu USclOel0Pu = 0 "Stator current overexcitation limitation initial output voltage in pu (base UNom)";
  parameter Types.VoltageModulePu USclUel0Pu = 0 "Stator current underexcitation limitation initial output voltage in pu (base UNom)";
  parameter Types.VoltageModulePu UUel0Pu = 0 "Underexcitation limitation initial output voltage in pu (base UNom)";

  //Initial parameter (calculated by initialization model)
  parameter Types.VoltageModulePu Vb0Pu "Initial available exciter field voltage in pu (base UNom)";

equation
  if PositionOel == 1 then
    sum1.u[2] = UOelPu;
    min1.u[2] = min1.u[1];
    sum2.u[3] = 0;
    min2.u[2] = min2.u[1];
  elseif PositionOel == 2 then
    sum1.u[2] = 0;
    min1.u[2] = UOelPu;
    sum2.u[3] = 0;
    min2.u[2] = min2.u[1];
  elseif PositionOel == 3 then
    sum1.u[2] = 0;
    min1.u[2] = min1.u[1];
    sum2.u[3] = UOelPu;
    min2.u[2] = min2.u[1];
  elseif PositionOel == 4 then
    sum1.u[2] = 0;
    min1.u[2] = min1.u[1];
    sum2.u[3] = 0;
    min2.u[2] = UOelPu;
  else
    sum1.u[2] = 0;
    min1.u[2] = min1.u[1];
    sum2.u[3] = 0;
    min2.u[2] = min2.u[1];
  end if;

  if PositionUel == 1 then
    sum1.u[3] = UUelPu;
    max1.u[2] = max1.u[1];
    sum2.u[4] = 0;
    max2.u[2] = max2.u[1];
  elseif PositionUel == 2 then
    sum1.u[3] = 0;
    max1.u[2] = UUelPu;
    sum2.u[4] = 0;
    max2.u[2] = max2.u[1];
  elseif PositionUel == 3 then
    sum1.u[3] = 0;
    max1.u[2] = max1.u[1];
    sum2.u[4] = UUelPu;
    max2.u[2] = max2.u[1];
  elseif PositionUel == 4 then
    sum1.u[3] = 0;
    max1.u[2] = max1.u[1];
    sum2.u[4] = 0;
    max2.u[2] = UUelPu;
  else
    sum1.u[3] = 0;
    max1.u[2] = max1.u[1];
    sum2.u[4] = 0;
    max2.u[2] = max2.u[1];
  end if;

  if PositionScl == 1 then
    sum1.u[4] = USclOelPu + USclUelPu;
    max1.u[3] = max1.u[1];
    min1.u[3] = min1.u[1];
    sum2.u[5] = 0;
    max2.u[3] = max2.u[1];
    min2.u[3] = min2.u[1];
  elseif PositionScl == 2 then
    sum1.u[4] = 0;
    max1.u[3] = USclUelPu;
    min1.u[3] = USclOelPu;
    sum2.u[5] = 0;
    max2.u[3] = max2.u[1];
    min2.u[3] = min2.u[1];
  elseif PositionScl == 3 then
    sum1.u[4] = 0;
    max1.u[3] = max1.u[1];
    min1.u[3] = min1.u[1];
    sum2.u[5] = USclOelPu + USclUelPu;
    max2.u[3] = max2.u[1];
    min2.u[3] = min2.u[1];
  elseif PositionScl == 4 then
    sum1.u[4] = 0;
    max1.u[3] = max1.u[1];
    min1.u[3] = min1.u[1];
    sum2.u[5] = 0;
    max2.u[3] = USclUelPu;
    min2.u[3] = USclOelPu;
  else
    sum1.u[4] = 0;
    max1.u[3] = max1.u[1];
    min1.u[3] = min1.u[1];
    sum2.u[5] = 0;
    max2.u[3] = max2.u[1];
    min2.u[3] = min2.u[1];
  end if;

  connect(const.y, min5.u1) annotation(
    Line(points = {{-59, 200}, {-40, 200}, {-40, 166}, {-22, 166}}, color = {0, 0, 127}));
  connect(product1.y, min5.u2) annotation(
    Line(points = {{-59, 140}, {-40, 140}, {-40, 154}, {-22, 154}}, color = {0, 0, 127}));
  connect(division.y, rectifierRegulationCharacteristic.u) annotation(
    Line(points = {{-159, 160}, {-142, 160}}, color = {0, 0, 127}));
  connect(UsPu, firstOrder.u) annotation(
    Line(points = {{-460, -80}, {-422, -80}}, color = {0, 0, 127}));
  connect(product.y, EfdPu) annotation(
    Line(points = {{381, 0}, {430, 0}}, color = {0, 0, 127}));
  connect(gain1.y, division.u1) annotation(
    Line(points = {{-219, 200}, {-200, 200}, {-200, 166}, {-182, 166}}, color = {0, 0, 127}));
  connect(IrPu, gain1.u) annotation(
    Line(points = {{-460, 200}, {-242, 200}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
  connect(limPI1.y, feedback.u1) annotation(
    Line(points = {{-39, -100}, {-8, -100}}, color = {0, 0, 127}));
  connect(rectifierRegulationCharacteristic.y, product1.u1) annotation(
    Line(points = {{-119, 160}, {-100, 160}, {-100, 146}, {-82, 146}}, color = {0, 0, 127}));
  connect(UsRefPu, add.u1) annotation(
    Line(points = {{-460, -40}, {-380, -40}, {-380, -54}, {-362, -54}}, color = {0, 0, 127}));
  connect(firstOrder.y, add.u2) annotation(
    Line(points = {{-399, -80}, {-380, -80}, {-380, -66}, {-362, -66}}, color = {0, 0, 127}));
  connect(limPI1.y, add2.u1) annotation(
    Line(points = {{-39, -100}, {-20, -100}, {-20, -74}, {38, -74}}, color = {0, 0, 127}));
  connect(feedback.y, add2.u2) annotation(
    Line(points = {{9, -100}, {20, -100}, {20, -86}, {38, -86}}, color = {0, 0, 127}));
  connect(add2.y, limiter.u) annotation(
    Line(points = {{61, -80}, {78, -80}}, color = {0, 0, 127}));
  connect(min4.y, limitedFirstOrder.u) annotation(
    Line(points = {{281, -20}, {298, -20}}, color = {0, 0, 127}));
  connect(product.y, firstOrder1.u) annotation(
    Line(points = {{381, 0}, {400, 0}, {400, -200}, {322, -200}}, color = {0, 0, 127}));
  connect(firstOrder1.y, feedback.u2) annotation(
    Line(points = {{299, -200}, {0, -200}, {0, -108}}, color = {0, 0, 127}));
  connect(utPu, potentialCircuit.uT) annotation(
    Line(points = {{-460, 160}, {-420, 160}, {-420, 144}, {-402, 144}}, color = {85, 170, 255}));
  connect(itPu, potentialCircuit.iT) annotation(
    Line(points = {{-460, 120}, {-420, 120}, {-420, 136}, {-402, 136}}, color = {85, 170, 255}));
  connect(min5.y, product.u1) annotation(
    Line(points = {{1, 160}, {340, 160}, {340, 6}, {358, 6}}, color = {0, 0, 127}));
  connect(limitedFirstOrder.y, product.u2) annotation(
    Line(points = {{321, -20}, {340, -20}, {340, -6}, {358, -6}}, color = {0, 0, 127}));
  connect(IrPu, feedback1.u2) annotation(
    Line(points = {{-460, 200}, {-280, 200}, {-280, 28}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
  connect(const1.y, feedback1.u1) annotation(
    Line(points = {{-319, 20}, {-288, 20}}, color = {0, 0, 127}));
  connect(feedback1.y, gain.u) annotation(
    Line(points = {{-271, 20}, {-242, 20}}, color = {0, 0, 127}));
  connect(max3.y, min4.u1) annotation(
    Line(points = {{-159, 40}, {240, 40}, {240, -14}, {258, -14}}, color = {0, 0, 127}));
  connect(max3.y, min3.u1) annotation(
    Line(points = {{-159, 40}, {-140, 40}, {-140, 26}, {-122, 26}}, color = {0, 0, 127}));
  connect(const3.y, min3.u2) annotation(
    Line(points = {{-159, 0}, {-140, 0}, {-140, 14}, {-122, 14}}, color = {0, 0, 127}));
  connect(min3.y, limPI1.limitMax) annotation(
    Line(points = {{-99, 20}, {-80, 20}, {-80, -94}, {-62, -94}}, color = {0, 0, 127}));
  connect(const4.y, limPI1.limitMin) annotation(
    Line(points = {{-99, -200}, {-80, -200}, {-80, -106}, {-62, -106}}, color = {0, 0, 127}));
  connect(UPssPu, sum2.u[2]) annotation(
    Line(points = {{-460, -200}, {-140, -200}, {-140, -100}, {-122, -100}}, color = {0, 0, 127}));
  connect(sum2.y, limPI1.u) annotation(
    Line(points = {{-99, -100}, {-63, -100}}, color = {0, 0, 127}));
  connect(add.y, sum1.u[1]) annotation(
    Line(points = {{-338, -60}, {-302, -60}}, color = {0, 0, 127}));
  connect(sum1.y, max1.u[1]) annotation(
    Line(points = {{-278, -60}, {-240, -60}}, color = {0, 0, 127}));
  connect(max1.yMax, min1.u[1]) annotation(
    Line(points = {{-218, -54}, {-180, -54}}, color = {0, 0, 127}));
  connect(min1.yMin, sum2.u[1]) annotation(
    Line(points = {{-158, -60}, {-140, -60}, {-140, -100}, {-122, -100}}, color = {0, 0, 127}));
  connect(limiter.y, max2.u[1]) annotation(
    Line(points = {{102, -80}, {140, -80}}, color = {0, 0, 127}));
  connect(max2.yMax, min2.u[1]) annotation(
    Line(points = {{162, -74}, {200, -74}}, color = {0, 0, 127}));
  connect(min2.yMax, min4.u2) annotation(
    Line(points = {{222, -68}, {240, -68}, {240, -26}, {258, -26}}, color = {0, 0, 127}));
  connect(const2.y, max3.u1) annotation(
    Line(points = {{-219, 60}, {-200, 60}, {-200, 46}, {-183, 46}}, color = {0, 0, 127}));
  connect(gain.y, max3.u2) annotation(
    Line(points = {{-219, 20}, {-200, 20}, {-200, 34}, {-183, 34}}, color = {0, 0, 127}));
  connect(potentialCircuit.vE, switch.u1) annotation(
    Line(points = {{-378, 140}, {-360, 140}, {-360, 108}, {-342, 108}}, color = {0, 0, 127}));
  connect(booleanConstant.y, switch.u2) annotation(
    Line(points = {{-378, 100}, {-342, 100}}, color = {255, 0, 255}));
  connect(const5.y, switch.u3) annotation(
    Line(points = {{-378, 60}, {-360, 60}, {-360, 92}, {-342, 92}}, color = {0, 0, 127}));
  connect(switch.y, division.u2) annotation(
    Line(points = {{-318, 100}, {-200, 100}, {-200, 154}, {-182, 154}}, color = {0, 0, 127}));
  connect(switch.y, product1.u2) annotation(
    Line(points = {{-318, 100}, {-100, 100}, {-100, 134}, {-82, 134}}, color = {0, 0, 127}));

  annotation(
    preferredView = "diagram",
    Diagram(coordinateSystem(extent = {{-440, -220}, {420, 220}})));
end St6c;
