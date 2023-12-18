within Dynawo.Electrical.Controls.WECC.BaseControls;

/*
* Copyright (c) 2021, RTE (http://www.rte-france.com)
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

model IqInjectionLogic "Reactive Current Injection Logic"

  //Parameters
  parameter Types.PerUnit IqFrzPu "Value at which reactive current injection (during a voltage-dip) is held for tHld seconds following a voltage dip if tHld > 0; in pu (base SNom, UNom) (typical: 0..Iqh1)";
  parameter Types.Time tHld "Time for which reactive current injection is held at some value following termination of the voltage-dip; if positive, then current is held at IqFrzPu, if negative then held at the value just prior to ending of the voltage-dip; in s (typical: -1..1, set to 0 to disable)";

  //Input variables
  Modelica.Blocks.Interfaces.RealInput iqVPu(start = 0) "Input for voltage-dependent reactive current injection in pu (base SNom, UNom)" annotation(
    Placement(visible = true, transformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.BooleanInput vDip(start = false) "True if voltage dip" annotation(
    Placement(visible = true, transformation(origin = {-120, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, 80}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));

  //Output variable
  Modelica.Blocks.Interfaces.RealOutput iqInjPu(start = 0) "Reactive current injection in pu (base SNom, UNom)" annotation(
    Placement(visible = true, transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

  Dynawo.NonElectrical.Blocks.MathBoolean.OffDelay holdAfterDip(tDelay = abs(tHld)) "True for |tHld| seconds after a voltage dip" annotation(
    Placement(visible = true, transformation(origin = {-76, 40}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));

equation
  if vDip or holdAfterDip.y and tHld < 0 then
    iqInjPu = iqVPu;
  elseif holdAfterDip.y and tHld > 0 then
    iqInjPu = IqFrzPu;
  else
    iqInjPu = 0;
  end if;

  connect(vDip, holdAfterDip.u) annotation(
    Line(points = {{-120, 40}, {-82, 40}}, color = {255, 0, 255}));

  annotation(
    preferredView = "text",
    Documentation(info = "<html><head></head><body><p>This block implements the behavior of the switch mechanic for reactive current injection, as specified in:<br><a href=\"https://www.wecc.org/Reliability/WECC-Second-Generation-Wind-Turbine-Models-012314.pdf\">https://www.wecc.org/Reliability/WECC-Second-Generation-Wind-Turbine-Models-012314.pdf</a></p><ul><li>Setting tHld to 0 results in abrupt ending of injection after the voltage dip has ended.</li><li>Setting tHld to a negative value continues the injection with the voltage-dependent injection for the absolute value of tHld seconds after the voltage dip has ended.</li><li>Setting tHld to a positive value continues the injection with a set constant (IqFrzPu) for the absolute value of tHld seconds after the voltage dip has ended.</li></ul></body><img src=\"modelica://Dynawo/Electrical/Controls/WECC/Images/IqInjStates.png\" alt=\"State transition diagram for the reec_a model\"></html>"),
    Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {-2, 74}, extent = {{-74, -38}, {82, -78}}, textString = "Reactive Current"), Text(origin = {1, -7}, extent = {{-63, 17}, {69, -21}}, textString = "Injection Logic"), Text(origin = {-129, 108}, extent = {{-19, 10}, {19, -10}}, textString = "vDip"), Text(origin = {-127, 28}, extent = {{-19, 10}, {19, -10}}, textString = "iqVPu"), Text(origin = {-127, 28}, extent = {{-19, 10}, {19, -10}}, textString = "iqVPu"), Text(origin = {121, 16}, extent = {{-19, 10}, {19, -10}}, textString = "iqInjPu")}));
end IqInjectionLogic;
