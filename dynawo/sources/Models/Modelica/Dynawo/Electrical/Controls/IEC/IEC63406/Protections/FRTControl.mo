within Dynawo.Electrical.Controls.IEC.IEC63406.Protections;

/*
* Copyright (c) 2024, RTE (http://www.rte-france.com)
* See AUTHORS.txt
* All rights reserved.
* This Source Code Form is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
* file, you can obtain one at http://mozilla.org/MPL/2.0/.
* SPDX-License-Identifier: MPL-2.0
*
* This file is part of Dynawo, an hybrid C++/Modelica open source suite of simulation tools for power systems.
*/

model FRTControl

  //Nominal parameter
  parameter Types.ApparentPowerModule SNom "Nominal converter apparent power in MVA";

  //General parameters
  parameter Types.PerUnit IPMaxPu "Maximum active current" annotation(
    Dialog(tab = "FRT"));
  parameter Types.PerUnit IQMaxPu "Maximum reactive current" annotation(
    Dialog(tab = "FRT"));
  parameter Types.PerUnit IPMinPu "Minimum active current" annotation(
    Dialog(tab = "FRT"));
  parameter Types.PerUnit IQMinPu "Minimum reactive current" annotation(
    Dialog(tab = "FRT"));
  parameter Types.PerUnit IMaxPu "Maximum current" annotation(
      Dialog(tab = "FRT"));
  parameter Boolean pqFRTFlag "Active/reactive control priority, 0/1" annotation(
      Dialog(tab = "FRT"));

  //LVRT and HVRT parameters
  parameter Types.PerUnit uLVRTPu "LVRT threshold value" annotation(
      Dialog(tab = "FRT"));
  parameter Types.PerUnit uHVRTPu "HVRT threshold value" annotation(
      Dialog(tab = "FRT"));
  parameter Boolean LVRTinPFlag "Active current flag during LVRT, 0/1" annotation(
      Dialog(tab = "FRT"));
  parameter Boolean LVRTinQFlag "Reactive current flag during LVRT, 0/1" annotation(
      Dialog(tab = "FRT"));
  parameter Boolean HVRTinPFlag "Active current flag during HVRT, 0/1" annotation(
      Dialog(tab = "FRT"));
  parameter Boolean HVRTinQFlag "Reactive current flag during HVRT, 0/1" annotation(
      Dialog(tab = "FRT"));

  parameter Real K1IpLV "Active current factor 1 during LVRT" annotation(
      Dialog(tab = "FRT"));
  parameter Real K2IpLV "Active current factor 2 during LVRT" annotation(
      Dialog(tab = "FRT"));
  parameter Real K1IqLV "Reactive current factor 1 during LVRT" annotation(
      Dialog(tab = "FRT"));
  parameter Real K2IqLV "Reactive current factor 2 during LVRT" annotation(
      Dialog(tab = "FRT"));
  parameter Real KpLVRT "Active power factor during LVRT" annotation(
      Dialog(tab = "FRT"));
  parameter Real KqLVRT "Reactive power factor during LVRT" annotation(
      Dialog(tab = "FRT"));
  parameter Types.PerUnit iPSetLVPu "Active current setting during LVRT" annotation(
      Dialog(tab = "FRT"));
  parameter Types.PerUnit iQSetLVPu "Reactive current setting during LVRT" annotation(
      Dialog(tab = "FRT"));
  parameter Types.PerUnit pSetLVPu "Active power setting during LVRT" annotation(
      Dialog(tab = "FRT"));
  parameter Types.PerUnit qSetLVPu "Reactive power setting during LVRT" annotation(
      Dialog(tab = "FRT"));

  parameter Real K1IpHV "Active current factor 1 during HVRT" annotation(
      Dialog(tab = "FRT"));
  parameter Real K2IpHV "Active current factor 2 during HVRT" annotation(
      Dialog(tab = "FRT"));
  parameter Real K1IqHV "Reactive current factor 1 during HVRT" annotation(
      Dialog(tab = "FRT"));
  parameter Real K2IqHV "Reactive current factor 2 during HVRT" annotation(
      Dialog(tab = "FRT"));
  parameter Real KpHVRT "Active power factor during HVRT" annotation(
      Dialog(tab = "FRT"));
  parameter Real KqHVRT "Reactive power factor during HVRT" annotation(
      Dialog(tab = "FRT"));
  parameter Types.PerUnit iPSetHVPu "Active current setting during HVRT" annotation(
      Dialog(tab = "FRT"));
  parameter Types.PerUnit iQSetHVPu "Reactive current setting during HVRT" annotation(
      Dialog(tab = "FRT"));
  parameter Types.PerUnit pSetHVPu "Active power setting during HVRT" annotation(
      Dialog(tab = "FRT"));
  parameter Types.PerUnit qSetHVPu "Reactive power setting during HVRT" annotation(
      Dialog(tab = "FRT"));


  //Input variables
  Modelica.Blocks.Interfaces.RealInput iPcmdPu(start = -P0Pu*SystemBase.SnRef/(SNom*U0Pu)) "Active current command" annotation(
      Placement(visible = true, transformation(origin = {-180, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-180, 120}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput iQcmdPu(start = Q0Pu*SystemBase.SnRef/(SNom*U0Pu)) "Reactive current command" annotation(
      Placement(visible = true, transformation(origin = {-180, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-180, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput uMeasPu(start = U0Pu) "Measured (and filtered) voltage component" annotation(
      Placement(visible = true, transformation(origin = {-180, 120}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-180, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput pMeasPu(start = -P0Pu * SystemBase.SnRef / SNom) annotation(
    Placement(visible = true, transformation(origin = {-180, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-180, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput qMeasPu(start = -Q0Pu * SystemBase.SnRef / SNom) annotation(
    Placement(visible = true, transformation(origin = {-180, -120}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-180, -120}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  //Output variables
  Modelica.Blocks.Interfaces.RealOutput ippPu(start = -P0Pu*SystemBase.SnRef/(SNom*U0Pu)) "Intermediary active current command in protection block" annotation(
      Placement(visible = true, transformation(origin = {180, 100}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {180, 100}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput iqqPu(start = Q0Pu*SystemBase.SnRef/(SNom*U0Pu)) "Intermediary reactive current command in protection block" annotation(
      Placement(visible = true, transformation(origin = {180, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {180, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput ipLVRTPu annotation(
      Placement(visible = true, transformation(origin = {170, 140}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {170, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput iqLVRTPu annotation(
      Placement(visible = true, transformation(origin = {170, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {170, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput ipHVRTPu annotation(
      Placement(visible = true, transformation(origin = {170, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {170, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput iqHVRTPu annotation(
      Placement(visible = true, transformation(origin = {170, -140}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {170, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

  //Initial parameters
  parameter Types.PerUnit IPMin0Pu annotation(
      Dialog(tab = "Operating point"));
  parameter Types.PerUnit IPMax0Pu annotation(
      Dialog(tab = "Operating point"));
  parameter Types.PerUnit IQMin0Pu annotation(
      Dialog(tab = "Operating point"));
  parameter Types.PerUnit IQMax0Pu annotation(
      Dialog(tab = "Operating point"));
  parameter Types.VoltageModulePu U0Pu "Initial voltage amplitude at grid terminal in pu (base UNom)" annotation(
    Dialog(group="Operating point"));
  parameter Types.ActivePowerPu P0Pu "Initial active power at grid terminal in pu (base SnRef) (receptor convention)" annotation(
    Dialog(tab = "Operating point"));
  parameter Types.ReactivePowerPu Q0Pu "Initial reactive power at grid terminal in pu (base SnRef) (receptor convention)" annotation(
    Dialog(tab = "Operating point"));

  Modelica.Blocks.Logical.Switch switch1 annotation(
      Placement(visible = true, transformation(origin = {-20, 140}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Blocks.Logical.Switch switch11 annotation(
      Placement(visible = true, transformation(origin = {-20, 20}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanExpression booleanExpression(y = LVRTinPFlag)  annotation(
      Placement(visible = true, transformation(origin = {-80, 140}, extent = {{-20, -10}, {20, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanExpression booleanExpression1(y = LVRTinQFlag) annotation(
      Placement(visible = true, transformation(origin = {-80, 20}, extent = {{-20, -10}, {20, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.VariableLimiter variableLimiter annotation(
      Placement(visible = true, transformation(origin = {140, 140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.VariableLimiter variableLimiter1 annotation(
      Placement(visible = true, transformation(origin = {140, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.VariableLimiter variableLimiter2 annotation(
      Placement(visible = true, transformation(origin = {140, -140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.Switch switch12 annotation(
      Placement(visible = true, transformation(origin = {-20, -140}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.VariableLimiter variableLimiter3 annotation(
      Placement(visible = true, transformation(origin = {140, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.Switch switch13 annotation(
      Placement(visible = true, transformation(origin = {-20, -20}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanExpression booleanExpression2(y = HVRTinQFlag) annotation(
      Placement(visible = true, transformation(origin = {-80, -140}, extent = {{-20, -10}, {20, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanExpression booleanExpression3(y = HVRTinPFlag) annotation(
      Placement(visible = true, transformation(origin = {-80, -20}, extent = {{-20, -10}, {20, 10}}, rotation = 0)));
  Dynawo.Electrical.Controls.IEC.IEC63406.Protections.AuxiliaryProtections.FRTCurrentBounds LVRTCurrentBounds(IMaxPu = IMaxPu, IPMaxPu = IPMaxPu, IPMinPu = IPMinPu, IQMaxPu = IQMaxPu, IQMinPu = IQMinPu, pqFRTFlag = pqFRTFlag)  annotation(
    Placement(visible = true, transformation(origin = {60, 80}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Dynawo.Electrical.Controls.IEC.IEC63406.Protections.AuxiliaryProtections.FRTCurrentBounds HVRTCurrentBounds1(IMaxPu = IMaxPu, IPMaxPu = IPMaxPu, IPMinPu = IPMinPu, IQMaxPu = IQMaxPu, IQMinPu = IQMinPu, pqFRTFlag = pqFRTFlag)  annotation(
    Placement(visible = true, transformation(origin = {60, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Dynawo.Electrical.Controls.IEC.IEC63406.Protections.AuxiliaryProtections.FRTCurrentCalculation LVRTCurrentCalculation(K1Ip = K1IpLV, K1Iq = K1IqLV, K2Ip = K2IpLV, K2Iq = K2IqLV, KpRT = KpLVRT, KqRT = KqLVRT, P0Pu = P0Pu, Q0Pu = Q0Pu, SNom = SNom, U0Pu = U0Pu, iPSetPu = iPSetLVPu, iQSetPu = iQSetLVPu, pSetPu = pSetLVPu, qSetPu = qSetLVPu, uHVRTPu = uHVRTPu, uLVRTPu = uLVRTPu, uRTPu = uLVRTPu)  annotation(
    Placement(visible = true, transformation(origin = {-80, 80}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Dynawo.Electrical.Controls.IEC.IEC63406.Protections.AuxiliaryProtections.FRTCurrentCalculation HVRTCurrentCalculation(K1Ip = K1IpHV, K1Iq = K1IqHV, K2Ip = K2IpHV, K2Iq = K2IqHV, KpRT = KpHVRT, KqRT = KqHVRT, P0Pu = P0Pu, Q0Pu = Q0Pu, SNom = SNom, U0Pu = U0Pu, iPSetPu = iPSetHVPu, iQSetPu = iQSetHVPu, pSetPu = pSetHVPu, qSetPu = qSetHVPu, uHVRTPu = uHVRTPu, uLVRTPu = uLVRTPu, uRTPu = uHVRTPu)  annotation(
    Placement(visible = true, transformation(origin = {-80, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
equation
  ippPu = if uMeasPu > uHVRTPu then ipHVRTPu else if uMeasPu < uLVRTPu then ipLVRTPu else iPcmdPu;
  iqqPu = if uMeasPu > uHVRTPu then iqHVRTPu else if uMeasPu < uLVRTPu then iqLVRTPu else iQcmdPu;
  connect(booleanExpression1.y, switch11.u2) annotation(
    Line(points = {{-58, 20}, {-32, 20}}, color = {255, 0, 255}));
  connect(booleanExpression.y, switch1.u2) annotation(
    Line(points = {{-58, 140}, {-31.8, 140}}, color = {255, 0, 255}));
  connect(switch11.y, variableLimiter1.u) annotation(
    Line(points = {{-9, 20}, {128, 20}}, color = {0, 0, 127}));
  connect(variableLimiter1.y, iqLVRTPu) annotation(
    Line(points = {{151, 20}, {170, 20}}, color = {0, 0, 127}));
  connect(switch1.y, variableLimiter.u) annotation(
    Line(points = {{-9, 140}, {128, 140}}, color = {0, 0, 127}));
  connect(variableLimiter.y, ipLVRTPu) annotation(
    Line(points = {{151, 140}, {170, 140}}, color = {0, 0, 127}));
  connect(switch12.y, variableLimiter2.u) annotation(
    Line(points = {{-9, -140}, {128, -140}}, color = {0, 0, 127}));
  connect(switch13.y, variableLimiter3.u) annotation(
    Line(points = {{-9, -20}, {128, -20}}, color = {0, 0, 127}));
  connect(variableLimiter3.y, ipHVRTPu) annotation(
    Line(points = {{151, -20}, {170, -20}}, color = {0, 0, 127}));
  connect(variableLimiter2.y, iqHVRTPu) annotation(
    Line(points = {{151, -140}, {170, -140}}, color = {0, 0, 127}));
  connect(booleanExpression2.y, switch12.u2) annotation(
    Line(points = {{-58, -140}, {-32, -140}}, color = {255, 0, 255}));
  connect(booleanExpression3.y, switch13.u2) annotation(
    Line(points = {{-58, -20}, {-32, -20}}, color = {255, 0, 255}));
  connect(uMeasPu, LVRTCurrentCalculation.uMeasPu) annotation(
    Line(points = {{-180, 120}, {-112, 120}, {-112, 96}, {-104, 96}}, color = {0, 0, 127}));
  connect(uMeasPu, HVRTCurrentCalculation.uMeasPu) annotation(
    Line(points = {{-180, 120}, {-112, 120}, {-112, -64}, {-104, -64}}, color = {0, 0, 127}));
  connect(iPcmdPu, LVRTCurrentCalculation.iPcmdPu) annotation(
    Line(points = {{-180, 60}, {-120, 60}, {-120, 88}, {-104, 88}}, color = {0, 0, 127}));
  connect(iPcmdPu, HVRTCurrentCalculation.iPcmdPu) annotation(
    Line(points = {{-180, 60}, {-120, 60}, {-120, -72}, {-104, -72}}, color = {0, 0, 127}));
  connect(iQcmdPu, LVRTCurrentCalculation.iQcmdPu) annotation(
    Line(points = {{-180, 0}, {-130, 0}, {-130, 80}, {-104, 80}}, color = {0, 0, 127}));
  connect(iQcmdPu, HVRTCurrentCalculation.iQcmdPu) annotation(
    Line(points = {{-180, 0}, {-130, 0}, {-130, -80}, {-104, -80}}, color = {0, 0, 127}));
  connect(pMeasPu, HVRTCurrentCalculation.pMeasPu) annotation(
    Line(points = {{-180, -60}, {-140, -60}, {-140, -88}, {-104, -88}}, color = {0, 0, 127}));
  connect(pMeasPu, LVRTCurrentCalculation.pMeasPu) annotation(
    Line(points = {{-180, -60}, {-140, -60}, {-140, 72}, {-104, 72}}, color = {0, 0, 127}));
  connect(qMeasPu, LVRTCurrentCalculation.qMeasPu) annotation(
    Line(points = {{-180, -120}, {-150, -120}, {-150, 64}, {-104, 64}}, color = {0, 0, 127}));
  connect(qMeasPu, HVRTCurrentCalculation.qMeasPu) annotation(
    Line(points = {{-180, -120}, {-150, -120}, {-150, -96}, {-104, -96}}, color = {0, 0, 127}));
  connect(LVRTCurrentCalculation.ipRTPu0, switch1.u3) annotation(
    Line(points = {{-58, 94}, {-48, 94}, {-48, 148}, {-32, 148}}, color = {0, 0, 127}));
  connect(LVRTCurrentCalculation.ipRTPu1, switch1.u1) annotation(
    Line(points = {{-58, 86}, {-40, 86}, {-40, 132}, {-32, 132}}, color = {0, 0, 127}));
  connect(LVRTCurrentCalculation.iqRTPu0, switch11.u3) annotation(
    Line(points = {{-58, 74}, {-40, 74}, {-40, 28}, {-32, 28}}, color = {0, 0, 127}));
  connect(LVRTCurrentCalculation.iqRTPu1, switch11.u1) annotation(
    Line(points = {{-58, 66}, {-48, 66}, {-48, 12}, {-32, 12}}, color = {0, 0, 127}));
  connect(HVRTCurrentCalculation.ipRTPu0, switch13.u3) annotation(
    Line(points = {{-58, -66}, {-50, -66}, {-50, -12}, {-32, -12}}, color = {0, 0, 127}));
  connect(HVRTCurrentCalculation.ipRTPu1, switch13.u1) annotation(
    Line(points = {{-58, -74}, {-40, -74}, {-40, -28}, {-32, -28}}, color = {0, 0, 127}));
  connect(HVRTCurrentCalculation.iqRTPu0, switch12.u3) annotation(
    Line(points = {{-58, -86}, {-40, -86}, {-40, -132}, {-32, -132}}, color = {0, 0, 127}));
  connect(HVRTCurrentCalculation.iqRTPu1, switch12.u1) annotation(
    Line(points = {{-58, -94}, {-50, -94}, {-50, -148}, {-32, -148}}, color = {0, 0, 127}));
  connect(switch13.y, HVRTCurrentBounds1.ipLVRTPrimPu) annotation(
    Line(points = {{-9, -20}, {0, -20}, {0, -72}, {36, -72}}, color = {0, 0, 127}));
  connect(switch12.y, HVRTCurrentBounds1.iqLVRTPrimPu) annotation(
    Line(points = {{-9, -140}, {0, -140}, {0, -88}, {36, -88}}, color = {0, 0, 127}));
  connect(switch11.y, LVRTCurrentBounds.iqLVRTPrimPu) annotation(
    Line(points = {{-9, 20}, {0, 20}, {0, 72}, {36, 72}}, color = {0, 0, 127}));
  connect(switch1.y, LVRTCurrentBounds.ipLVRTPrimPu) annotation(
    Line(points = {{-9, 140}, {0, 140}, {0, 88}, {36, 88}}, color = {0, 0, 127}));
  connect(LVRTCurrentBounds.ipMaxPu, variableLimiter.limit1) annotation(
    Line(points = {{82, 94}, {110, 94}, {110, 148}, {128, 148}}, color = {0, 0, 127}));
  connect(LVRTCurrentBounds.ipMinPu, variableLimiter.limit2) annotation(
    Line(points = {{82, 86}, {120, 86}, {120, 132}, {128, 132}}, color = {0, 0, 127}));
  connect(LVRTCurrentBounds.iqMaxPu, variableLimiter1.limit1) annotation(
    Line(points = {{82, 74}, {120, 74}, {120, 28}, {128, 28}}, color = {0, 0, 127}));
  connect(LVRTCurrentBounds.iqMinPu, variableLimiter1.limit2) annotation(
    Line(points = {{82, 66}, {110, 66}, {110, 12}, {128, 12}}, color = {0, 0, 127}));
  connect(HVRTCurrentBounds1.ipMaxPu, variableLimiter3.limit1) annotation(
    Line(points = {{82, -66}, {110, -66}, {110, -12}, {128, -12}}, color = {0, 0, 127}));
  connect(HVRTCurrentBounds1.ipMinPu, variableLimiter3.limit2) annotation(
    Line(points = {{82, -74}, {120, -74}, {120, -28}, {128, -28}}, color = {0, 0, 127}));
  connect(HVRTCurrentBounds1.iqMaxPu, variableLimiter2.limit1) annotation(
    Line(points = {{82, -86}, {120, -86}, {120, -132}, {128, -132}}, color = {0, 0, 127}));
  connect(HVRTCurrentBounds1.iqMinPu, variableLimiter2.limit2) annotation(
    Line(points = {{82, -94}, {110, -94}, {110, -148}, {128, -148}}, color = {0, 0, 127}));
  annotation(
    Icon(graphics = {Rectangle(extent = {{-160, 160}, {160, -160}}), Text(extent = {{-160, 160}, {160, -160}}, textString = "FRT
Control")}, coordinateSystem(extent = {{-160, -160}, {160, 160}})),
    Diagram(coordinateSystem(extent = {{-160, -160}, {160, 160}})));
end FRTControl;
