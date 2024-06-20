within Dynawo.Electrical.Controls.IEC.BaseClasses;

/*
* Copyright (c) 2022, RTE (http://www.rte-france.com)
* See AUTHORS.txt
* All rights reserved.
* This Source Code Form is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
* file, you can obtain one at http://mozilla.org/MPL/2.0/.
* SPDX-License-Identifier: MPL-2.0
*
* This file is part of Dynawo, an hybrid C++/Modelica open source suite of simulation tools for power systems.
*/

model BaseMeasurements "Measurement module for wind turbine controls (IEC N°61400-27-1)"

  //Nominal parameters
  parameter Types.ApparentPowerModule SNom "Nominal converter apparent power in MVA";
  parameter Types.Time tS "Integration time step in s";

  //Measurement parameters
  parameter Types.Time tIFilt "Filter time constant for current measurement in s" annotation(
    Dialog(tab = "Measurement"));
  parameter Types.Time tPFilt "Filter time constant for active power measurement in s" annotation(
    Dialog(tab = "Measurement"));
  parameter Types.Time tQFilt "Filter time constant for reactive power measurement in s" annotation(
    Dialog(tab = "Measurement"));
  parameter Types.Time tUFilt "Filter time constant for voltage measurement in s" annotation(
    Dialog(tab = "Measurement"));

  //Input variables
  Modelica.ComplexBlocks.Interfaces.ComplexInput iPu(re(start = -i0Pu.re * SystemBase.SnRef / SNom), im(start = -i0Pu.im * SystemBase.SnRef / SNom)) "Complex current at grid terminal in pu (base UNom, SNom) (generator convention)" annotation(
    Placement(visible = true, transformation(origin = {-160, 80}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.ComplexBlocks.Interfaces.ComplexInput uPu(re(start = u0Pu.re), im(start = u0Pu.im)) "Complex voltage at grid terminal in pu (base UNom)" annotation(
    Placement(visible = true, transformation(origin = {-160, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-110, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

  //Output variables
  Modelica.Blocks.Interfaces.RealOutput IFiltPu(start = ComplexMath.'abs'(i0Pu) * SystemBase.SnRef / SNom) "Filtered current module at grid terminal in pu (base UNom, SNom)" annotation(
    Placement(visible = true, transformation(origin = {150, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput PFiltPu(start = -P0Pu * SystemBase.SnRef / SNom) "Filtered active power at grid terminal in pu (base SNom) (generator convention)" annotation(
    Placement(visible = true, transformation(origin = {150, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput QFiltPu(start = -Q0Pu * SystemBase.SnRef / SNom) "Filtered reactive power at grid terminal in pu (base SNom) (generator convention)" annotation(
    Placement(visible = true, transformation(origin = {150, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput UFiltPu(start = U0Pu) "Filtered voltage amplitude at grid terminal in pu (base UNom)" annotation(
    Placement(visible = true, transformation(origin = {150, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

  Modelica.Blocks.Continuous.FirstOrder firstOrder(T = tPFilt, y_start = -P0Pu * SystemBase.SnRef / SNom) annotation(
    Placement(visible = true, transformation(origin = {90, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.FirstOrder firstOrder1(T = tQFilt, y_start = -Q0Pu * SystemBase.SnRef / SNom) annotation(
    Placement(visible = true, transformation(origin = {90, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.FirstOrder firstOrder2(T = tIFilt, y_start = ComplexMath.'abs'(i0Pu) * SystemBase.SnRef / SNom) annotation(
    Placement(visible = true, transformation(origin = {90, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.FirstOrder firstOrder3(T = tUFilt, y_start = U0Pu) annotation(
    Placement(visible = true, transformation(origin = {90, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.ComplexBlocks.ComplexMath.Product product(useConjugateInput2 = true) annotation(
    Placement(visible = true, transformation(origin = {-70, 100}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.ComplexBlocks.ComplexMath.ComplexToReal complexToReal annotation(
    Placement(visible = true, transformation(origin = {-10, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Dynawo.NonElectrical.Blocks.Complex.ComplexToPolar complexToPolar annotation(
    Placement(visible = true, transformation(origin = {-10, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Dynawo.NonElectrical.Blocks.Complex.ComplexToPolar complexToPolar1 annotation(
    Placement(visible = true, transformation(origin = {-70, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

  //Initial parameters
  parameter Types.ComplexCurrentPu i0Pu "Initial complex current at grid terminal in pu (base UNom, SnRef) (receptor convention)" annotation(
    Dialog(group = "Initialization"));
  parameter Types.ActivePowerPu P0Pu "Initial active power at grid terminal in pu (base SnRef) (receptor convention)" annotation(
    Dialog(tab = "Operating point"));
  parameter Types.ReactivePowerPu Q0Pu "Initial reactive power at grid terminal in pu (base SnRef) (receptor convention)" annotation(
    Dialog(tab = "Operating point"));
  parameter Types.VoltageModulePu U0Pu "Initial voltage amplitude at grid terminal in pu (base UNom)" annotation(
    Dialog(tab = "Operating point"));
  parameter Types.ComplexVoltagePu u0Pu "Initial complex voltage at grid terminal in pu (base UNom)" annotation(
    Dialog(group = "Initialization"));
  parameter Types.Angle UPhase0 "Initial voltage angle at grid terminal in rad" annotation(
    Dialog(tab = "Operating point"));


equation
  connect(iPu, product.u2) annotation(
    Line(points = {{-160, 80}, {-120, 80}, {-120, 106}, {-82, 106}}, color = {85, 170, 255}));
  connect(uPu, product.u1) annotation(
    Line(points = {{-160, 0}, {-100, 0}, {-100, 94}, {-82, 94}}, color = {85, 170, 255}));
  connect(product.y, complexToReal.u) annotation(
    Line(points = {{-58, 100}, {-22, 100}}, color = {85, 170, 255}));
  connect(complexToReal.re, firstOrder.u) annotation(
    Line(points = {{2, 106}, {40, 106}, {40, 120}, {78, 120}}, color = {0, 0, 127}));
  connect(complexToReal.im, firstOrder1.u) annotation(
    Line(points = {{2, 94}, {40, 94}, {40, 80}, {78, 80}}, color = {0, 0, 127}));
  connect(firstOrder.y, PFiltPu) annotation(
    Line(points = {{102, 120}, {150, 120}}, color = {0, 0, 127}));
  connect(firstOrder1.y, QFiltPu) annotation(
    Line(points = {{102, 80}, {150, 80}}, color = {0, 0, 127}));
  connect(iPu, complexToPolar.u) annotation(
    Line(points = {{-160, 80}, {-120, 80}, {-120, 40}, {-22, 40}}, color = {85, 170, 255}));
  connect(uPu, complexToPolar1.u) annotation(
    Line(points = {{-160, 0}, {-82, 0}}, color = {85, 170, 255}));
  connect(complexToPolar1.len, firstOrder3.u) annotation(
    Line(points = {{-58, 6}, {40, 6}, {40, -40}, {78, -40}}, color = {0, 0, 127}));
  connect(firstOrder3.y, UFiltPu) annotation(
    Line(points = {{102, -40}, {150, -40}}, color = {0, 0, 127}));
  connect(complexToPolar.len, firstOrder2.u) annotation(
    Line(points = {{2, 46}, {40, 46}, {40, 40}, {78, 40}}, color = {0, 0, 127}));
  connect(firstOrder2.y, IFiltPu) annotation(
    Line(points = {{102, 40}, {150, 40}}, color = {0, 0, 127}));

  annotation(
    preferredView = "diagram",
    Diagram(coordinateSystem(extent = {{-140, -140}, {140, 140}})),
    Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {100, 100}}), Text(origin = {16, -23}, extent = {{-108, -24}, {76, 10}}, textString = "Module"), Text(origin = {8, 35}, extent = {{-100, -30}, {88, 20}}, textString = "Measurement")}));
end BaseMeasurements;
