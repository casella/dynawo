within Dynawo.NonElectrical.Blocks.NonLinear.BaseClasses;

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

partial block BaseRSFlipFlop "Base block of RS flip flop"

  import Modelica;
  import Modelica.Blocks.Interfaces;

  extends Modelica.Blocks.Icons.PartialBooleanBlock;

  parameter Boolean yStart = false "Start value of y";

  //Input variables
  Interfaces.BooleanInput R annotation (
      Placement(transformation(extent={{-140,-80},{-100,-40}})));
  Interfaces.BooleanInput S annotation (
      Placement(transformation(extent={{-140,40},{-100,80}})));

  //Output variable
  Interfaces.BooleanOutput y(start = yStart) annotation (
      Placement(transformation(extent={{100,-10},{120,10}})));

  annotation (
  Icon(graphics={
      Text(
        extent={{-60,-30},{-20,-90}},
        textString="R"),
      Text(extent = {{-60, 90}, {-20, 30}}, textString = "S"),
      Ellipse(lineColor = {235, 235, 235}, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-73, 54}, {-87, 68}}, endAngle = 360),
      Ellipse(lineColor = {235, 235, 235}, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{-71, -52}, {-85, -66}}, endAngle = 360),
      Ellipse(lineColor = {235, 235, 235}, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid, extent = {{71, 7}, {85, -7}}, endAngle = 360)}),
    Documentation(info= "<html><head></head><body><p><br></p>
</body></html>"));
end BaseRSFlipFlop;
