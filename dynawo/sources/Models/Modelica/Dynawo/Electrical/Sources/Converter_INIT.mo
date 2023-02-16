within Dynawo.Electrical.Sources;

/*
* Copyright (c) 2015-2019, RTE (http://www.rte-france.com)
* See AUTHORS.txt
* All rights reserved.
* This Source Code Form is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
* file, you can obtain one at http://mozilla.org/MPL/2.0/.
* SPDX-License-Identifier: MPL-2.0
*
* This file is part of Dynawo, an hybrid C++/Modelica open source time domain simulation tool for power systems.
*/

model Converter_INIT "Initialization model for converter model for grid forming and grid following applications"
  import Modelica;
  import Modelica.ComplexMath;
  import Dynawo.Types;
  import Dynawo.Electrical.SystemBase;
  extends AdditionalIcons.Init;

  parameter Types.PerUnit RFilter "Filter resistance in pu (base UNom, SNom)";
  parameter Types.PerUnit LFilter "Filter inductance in pu (base UNom, SNom)";
  parameter Types.PerUnit CFilter "Filter capacitance in pu (base UNom, SNom)";
  parameter Types.PerUnit RTransformer "Transformer resistance in pu (base UNom, SNom)";
  parameter Types.PerUnit LTransformer "Transformer inductance in pu (base UNom, SNom)";
  parameter Types.ApparentPowerModule SNom "Apparent power module reference for the converter";

  parameter Types.VoltageModulePu U0Pu "Start value of voltage amplitude at terminal in pu (base UNom)";
  parameter Types.Angle UPhase0 "Start value of voltage angle at terminal in rad";
  parameter Types.ActivePowerPu P0Pu "Start value of active power in pu (base SnRef) (receptor convention)";
  parameter Types.ReactivePowerPu Q0Pu "Start value of reactive power in pu (base SnRef) (receptor convention)";

  Types.ComplexPerUnit i0Pu "Start value of the complex current at terminal in pu (base UNom, SnRef) (receptor convention)";
  Types.ComplexPerUnit u0Pu "Start value of the complex voltage at terminal in pu (base UNom)";
  Types.PerUnit UdPcc0Pu "Start value of the d-axis voltage at the PCC in pu (base UNom)";
  Types.PerUnit UqPcc0Pu "Start value of the q-axis voltage at the PCC in pu (base UNom)";
  Modelica.Blocks.Interfaces.RealOutput IdcSource0Pu "Start value of DC current in pu (base UdcNom, SNom)";
  Modelica.Blocks.Interfaces.RealInput UdcSource0Pu "Start value of DC voltage in pu (base UdcNom)";
  Modelica.Blocks.Interfaces.RealInput UdConv0Pu "Start value of d-axis modulation voltage in pu (base UNom)";
  Modelica.Blocks.Interfaces.RealInput UqConv0Pu "Start value of q-axis modulation voltage in pu (base UNom)";
  Modelica.Blocks.Interfaces.RealInput IdConv0Pu "Start value of d-axis current in the converter in pu (base UNom, SNom) (generator convention)";
  Modelica.Blocks.Interfaces.RealOutput IqConv0Pu "Start value of q-axis current in the converter in pu (base UNom, SNom) (generator convention)";
  Modelica.Blocks.Interfaces.RealOutput IdPcc0Pu "Start value of d-axis current in the grid in pu (base UNom, SNom) (generator convention)";
  Modelica.Blocks.Interfaces.RealInput IqPcc0Pu "Start value of q-axis current in the grid in pu (base UNom, SNom) (generator convention)";
  Modelica.Blocks.Interfaces.RealInput UdFilter0Pu "Start value of d-axis voltage at the converter's capacitor in pu (base UNom)";
  Modelica.Blocks.Interfaces.RealInput UqFilter0Pu "Start value of q-axis voltage at the converter's capacitor in pu (base UNom)";
  Modelica.Blocks.Interfaces.RealOutput Theta0 "Start value of phase shift between the converter's rotating frame and the grid rotating frame in rad";
  Modelica.Blocks.Interfaces.RealOutput PFilter0Pu "Start value of active power generated at the converter's capacitor in pu (base SNom) (generator convention)";
  Modelica.Blocks.Interfaces.RealOutput QFilter0Pu "Start value of reactive power generated at the converter's capacitor in pu (base SNom) (generator convention)";
  Types.PerUnit IConv0Pu "Start value of module of the current injected by the converter in pu (base UNom, SNom)";
  Types.PerUnit PGen0Pu "Start value of active power generated by the converter at the PCC in pu (base UNom, SnRef) (generator convention)";
  Types.PerUnit QGen0Pu "Start value of reactive power generated by the converter at the PCC in pu (base UNom, SnRef) (generator convention)";

equation
  u0Pu = ComplexMath.fromPolar(U0Pu, UPhase0);
  Complex(P0Pu, Q0Pu) = u0Pu * ComplexMath.conj(i0Pu);

  /* DQ reference frame change from network reference to converter reference and pu base change */
  [UdPcc0Pu; UqPcc0Pu] = [cos(Theta0), sin(Theta0); -sin(Theta0), cos(Theta0)] * [u0Pu.re; u0Pu.im];
  [IdPcc0Pu; IqPcc0Pu] = - [cos(Theta0), sin(Theta0); -sin(Theta0), cos(Theta0)] * [i0Pu.re; i0Pu.im] * SystemBase.SnRef / SNom;

  /* RL Transformer */
  0 = UdFilter0Pu - RTransformer * IdPcc0Pu + SystemBase.omegaRef0Pu * LTransformer * IqPcc0Pu - UdPcc0Pu;
  0 = UqFilter0Pu - RTransformer * IqPcc0Pu - SystemBase.omegaRef0Pu * LTransformer * IdPcc0Pu - UqPcc0Pu;

  /* RLC Filter */
  0 = UdConv0Pu - RFilter * IdConv0Pu + SystemBase.omegaRef0Pu * LFilter * IqConv0Pu - UdFilter0Pu;
  0 = UqConv0Pu - RFilter * IqConv0Pu - SystemBase.omegaRef0Pu * LFilter * IdConv0Pu - UqFilter0Pu;
  0 = IdConv0Pu + SystemBase.omegaRef0Pu * CFilter * UqFilter0Pu - IdPcc0Pu;
  0 = IqConv0Pu - SystemBase.omegaRef0Pu * CFilter * UdFilter0Pu - IqPcc0Pu;
  IConv0Pu = sqrt(IdConv0Pu * IdConv0Pu + IqConv0Pu * IqConv0Pu);

  /* Power Conservation */
  UdConv0Pu * IdConv0Pu + UqConv0Pu * IqConv0Pu = UdcSource0Pu * IdcSource0Pu;

  /* Power Calculation */
  PGen0Pu = (UdPcc0Pu * IdPcc0Pu + UqPcc0Pu * IqPcc0Pu) * SNom / SystemBase.SnRef;
  QGen0Pu = (UqPcc0Pu * IdPcc0Pu - UdPcc0Pu * IqPcc0Pu) * SNom / SystemBase.SnRef;
  PFilter0Pu = UdFilter0Pu * IdPcc0Pu + UqFilter0Pu * IqPcc0Pu;
  QFilter0Pu = UqFilter0Pu * IdPcc0Pu - UdFilter0Pu * IqPcc0Pu;

  annotation(preferredView = "text");
end Converter_INIT;
