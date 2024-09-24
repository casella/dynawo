within Dynawo.Electrical.Controls.IEC.IEC63406.Protections.AuxiliaryProtections;

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

model FRTCurrentCalculation

  //Nominal parameter
  parameter Types.ApparentPowerModule SNom "Nominal converter apparent power in MVA";

  //Parameters
  parameter Types.PerUnit uLVRTPu "LVRT threshold value" annotation(
    Dialog(tab = "FRT"));
  parameter Types.PerUnit uHVRTPu "HVRT threshold value" annotation(
    Dialog(tab = "FRT"));
  parameter Real K1Ip "Active current factor 1 during LVRT or HVRT" annotation(
    Dialog(tab = "FRT"));
  parameter Real K2Ip "Active current factor 2 during LVRT or HVRT" annotation(
    Dialog(tab = "FRT"));
  parameter Real K1Iq "Reactive current factor 1 during LVRT or HVRT" annotation(
    Dialog(tab = "FRT"));
  parameter Real K2Iq "Reactive current factor 2 during LVRT or HVRT" annotation(
    Dialog(tab = "FRT"));
  parameter Real KpRT "Active power factor during LVRT or HVRT" annotation(
    Dialog(tab = "FRT"));
  parameter Real KqRT "Reactive power factor during LVRT or HVRT" annotation(
    Dialog(tab = "FRT"));
  parameter Types.PerUnit iPSetPu "Active current setting during LVRT or HVRT" annotation(
    Dialog(tab = "FRT"));
  parameter Types.PerUnit iQSetPu "Reactive current setting during LVRT or HVRT" annotation(
    Dialog(tab = "FRT"));
  parameter Types.PerUnit pSetPu "Active power setting during LVRT or HVRT" annotation(
    Dialog(tab = "FRT"));
  parameter Types.PerUnit qSetPu "Reactive power setting during LVRT or HVRT" annotation(
    Dialog(tab = "FRT"));
  parameter Types.PerUnit uRTPu "LVRT or HVRT threshold value" annotation(
    Dialog(tab = "FRT"));
  Types.PerUnit pPreFaultPu(start = 0);
  Types.PerUnit qPreFaultPu(start = 0);

  //Input variables
  Modelica.Blocks.Interfaces.RealInput iPcmdPu(start = -P0Pu * SystemBase.SnRef / (SNom * U0Pu)) "Active current command" annotation(
    Placement(visible = true, transformation(origin = {-120, 80}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput iQcmdPu(start = Q0Pu * SystemBase.SnRef / (SNom * U0Pu)) "Reactive current command" annotation(
    Placement(visible = true, transformation(origin = {-120, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput uMeasPu(start = U0Pu) "Filtered voltage amplitude at grid terminal in pu (base UNom)" annotation(
    Placement(visible = true, transformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, 80}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput pMeasPu(start = -P0Pu * SystemBase.SnRef / SNom) "Filtered active power at grid terminal in pu (base SNom) (generator convention)" annotation(
    Placement(visible = true, transformation(origin = {-120, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput qMeasPu(start = -Q0Pu * SystemBase.SnRef / SNom) "Filtered reactive power at grid terminal in pu (base SNom) (generator convention)" annotation(
    Placement(visible = true, transformation(origin = {-120, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  //Output variables
  Modelica.Blocks.Interfaces.RealOutput ipRTPu0 annotation(
    Placement(visible = true, transformation(origin = {110, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput ipRTPu1 annotation(
    Placement(visible = true, transformation(origin = {110, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput iqRTPu0 annotation(
    Placement(visible = true, transformation(origin = {110, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput iqRTPu1 annotation(
    Placement(visible = true, transformation(origin = {110, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

  //Initial parameters
  parameter Types.VoltageModulePu U0Pu "Initial voltage amplitude at grid terminal in pu (base UNom)" annotation(
    Dialog(group = "Operating point"));
  parameter Types.ActivePowerPu P0Pu "Initial active power at grid terminal in pu (base SnRef) (receptor convention)" annotation(
    Dialog(tab = "Operating point"));
  parameter Types.ReactivePowerPu Q0Pu "Initial reactive power at grid terminal in pu (base SnRef) (receptor convention)" annotation(
    Dialog(tab = "Operating point"));

equation
  //PreFault powers for the fault control loop
  when uMeasPu < uLVRTPu then
    pPreFaultPu = pre(pMeasPu);
    qPreFaultPu = pre(qMeasPu);
  elsewhen uMeasPu > uHVRTPu then
    pPreFaultPu = pre(pMeasPu);
    qPreFaultPu = pre(qMeasPu);
  end when;

  ipRTPu0 = K1Ip * uMeasPu + K2Ip * iPcmdPu + iPSetPu;
  ipRTPu1 = (KpRT * pPreFaultPu + pSetPu) / uMeasPu;
  iqRTPu0 = K1Iq * (uRTPu - uMeasPu) + K2Iq * iQcmdPu + iQSetPu;
  iqRTPu1 = (KqRT * qPreFaultPu + qSetPu) / uMeasPu;

  annotation(
    Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}), Text(extent = {{-100, 100}, {100, -100}}, textString = "FRT
Current
Calculat°")}));
end FRTCurrentCalculation;
