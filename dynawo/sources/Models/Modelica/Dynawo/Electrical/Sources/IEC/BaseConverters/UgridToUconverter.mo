within Dynawo.Electrical.Sources.IEC.BaseConverters;

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

model UgridToUconverter

  //Parameters
  parameter Types.PerUnit XesPu "Coupling inductance in the filter or the transformer at the grid side" annotation(
    Dialog(tab = "Electrical"));
  parameter Types.PerUnit ResPu "Coupling resistance in the filter or the transformer at the grid side" annotation(
    Dialog(tab = "Electrical"));

  //Input variables
  Modelica.Blocks.Interfaces.RealInput idPu annotation(
    Placement(visible = true, transformation(origin = {-120, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput iqPu annotation(
    Placement(visible = true, transformation(origin = {-120, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput ugqPu annotation(
    Placement(visible = true, transformation(origin = {-40, -120}, extent = {{-20, -20}, {20, 20}}, rotation = 90), iconTransformation(origin = {-40, -120}, extent = {{-20, -20}, {20, 20}}, rotation = 90)));
  Modelica.Blocks.Interfaces.RealInput ugdPu annotation(
    Placement(visible = true, transformation(origin = {40, -120}, extent = {{-20, -20}, {20, 20}}, rotation = 90), iconTransformation(origin = {40, -120}, extent = {{-20, -20}, {20, 20}}, rotation = 90)));

  //Output variable
  Modelica.Blocks.Interfaces.RealOutput uedPu annotation(
    Placement(visible = true, transformation(origin = {110, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput ueqPu annotation(
    Placement(visible = true, transformation(origin = {110, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

equation
  uedPu = ugdPu + idPu * ResPu - iqPu * XesPu;
  ueqPu = ugqPu + iqPu * ResPu + idPu * XesPu;

annotation(
    Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}), Text(origin = {1, 0}, extent = {{-97, 96}, {97, -96}}, textString = "Ugrid\nTo\nUconverter")}));
end UgridToUconverter;
