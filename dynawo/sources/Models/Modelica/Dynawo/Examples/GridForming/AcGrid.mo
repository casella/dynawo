within Dynawo.Examples.GridForming;

model AcGrid "AC Grid from IEE explorer paper from Carmen C."
  
//  final parameter Real K_FH = 1 + FH;
  parameter Real SNom =SystemBase.SnRef;
  parameter Real U0pu;
  parameter Real UPhase0;
  parameter Real Upu(start=U0pu)  ;
  parameter Real UPhase(start=UPhase0) ;
  parameter Real omegaRefPu=SystemBase.omegaRef0Pu;
  

  parameter Real R = 0.05 "Governor Droop";
  parameter Real FH = 0.3 "Fraction of the total power generated by the HP turbine";
  parameter Real H = 4 "Inertie constant";
  parameter Real TR = 8 "Reheat time constant, seconds";
  parameter Real Km = 0.95 "Mechanical power gain factor";
  parameter Real D = 1 "Damping Factor";

  
  Modelica.Blocks.Continuous.FirstOrder firstOrder(T = TR, k = 1 - FH) annotation(
    Placement(visible = true, transformation(origin = {-38, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain FH_(k = FH) annotation(
    Placement(visible = true, transformation(origin = {-36, 68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add add annotation(
    Placement(visible = true, transformation(origin = {10, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add add1(k2 = -1) annotation(
    Placement(visible = true, transformation(origin = {-104, 68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain(k = 1 / R) annotation(
    Placement(visible = true, transformation(origin = {-27, -35}, extent = {{15, -15}, {-15, 15}}, rotation = 0)));
  Modelica.Blocks.Continuous.Integrator integrator(k = 1 / (2 * H)) annotation(
    Placement(visible = true, transformation(origin = {162, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain1(k = Km) annotation(
    Placement(visible = true, transformation(origin = {52, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain2(k = D) annotation(
    Placement(visible = true, transformation(origin = {138, 8}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add3 add3(k2 = -1, k3 = -1) annotation(
    Placement(visible = true, transformation(origin = {106, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.Integrator integrator1(k = 1) annotation(
    Placement(visible = true, transformation(origin = {396, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Step Psetpoint(height = 0, offset = 0)  annotation(
    Placement(visible = true, transformation(origin = {-182, 84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add add2 annotation(
    Placement(visible = true, transformation(origin = {274, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add add4(k2 = -1)  annotation(
    Placement(visible = true, transformation(origin = {340, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Dynawo.Examples.GridForming.PhaseurGrid phaseurGrid(SNom = SNom, UPhase = UPhase, UPhase0 = UPhase0, UPu = Upu, U0Pu = U0pu) annotation(
    Placement(visible = true, transformation(origin = {540, -18}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
  Dynawo.Connectors.ACPower aCPower annotation(
    Placement(visible = true, transformation(origin = {749, -17}, extent = {{-23, -23}, {23, 23}}, rotation = 0), iconTransformation(origin = {120, 74}, extent = {{-20, -20}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput omegaPu annotation(
    Placement(visible = true, transformation(origin = {414, 40}, extent = {{-15, -15}, {15, 15}}, rotation = 0), iconTransformation(extent = {{99, -73}, {129, -43}}, rotation = 0)));
  Modelica.Blocks.Sources.Step Wref(height = omegaRefPu) annotation(
    Placement(visible = true, transformation(origin = {222, -8}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain3(k = 1) annotation(
    Placement(visible = true, transformation(origin = {133, -35}, extent = {{15, -15}, {-15, 15}}, rotation = 0)));
equation
  connect(firstOrder.y, add.u2) annotation(
    Line(points = {{-26, 24}, {-14, 24}, {-14, 46}, {-2, 46}}, color = {0, 0, 127}));
  connect(add1.y, FH_.u) annotation(
    Line(points = {{-92, 68}, {-48, 68}}, color = {0, 0, 127}));
  connect(firstOrder.u, add1.y) annotation(
    Line(points = {{-50, 24}, {-80, 24}, {-80, 58}, {-92, 58}, {-92, 68}}, color = {0, 0, 127}));
  connect(add.y, gain1.u) annotation(
    Line(points = {{22, 52}, {40, 52}}, color = {0, 0, 127}));
  connect(gain1.y, add3.u1) annotation(
    Line(points = {{64, 52}, {74, 52}, {74, 58}, {94, 58}}, color = {0, 0, 127}));
  connect(add3.y, integrator.u) annotation(
    Line(points = {{118, 50}, {150, 50}}, color = {0, 0, 127}));
  connect(gain2.y, add3.u3) annotation(
    Line(points = {{128, 8}, {82, 8}, {82, 42}, {94, 42}}, color = {0, 0, 127}));
  connect(gain2.u, integrator.y) annotation(
    Line(points = {{150, 8}, {192, 8}, {192, 50}, {174, 50}}, color = {0, 0, 127}));
  connect(add1.y, FH_.u) annotation(
    Line(points = {{-92, 68}, {-48, 68}}, color = {0, 0, 127}));
  connect(gain.y, add1.u2) annotation(
    Line(points = {{-43.5, -35}, {-156, -35}, {-156, 62}, {-116, 62}}, color = {0, 0, 127}));
  connect(FH_.y, add.u1) annotation(
    Line(points = {{-24, 68}, {-18, 68}, {-18, 58}, {-2, 58}}, color = {0, 0, 127}));
  connect(Psetpoint.y, add1.u1) annotation(
    Line(points = {{-170, 84}, {-160, 84}, {-160, 74}, {-116, 74}}, color = {0, 0, 127}));
  connect(integrator.y, add2.u1) annotation(
    Line(points = {{174, 50}, {228, 50}, {228, 46}, {262, 46}}, color = {0, 0, 127}));
  connect(add4.u1, add2.y) annotation(
    Line(points = {{328, -38}, {284, -38}, {284, -2}, {314, -2}, {314, 40}, {286, 40}}, color = {0, 0, 127}));
  connect(add4.y, integrator1.u) annotation(
    Line(points = {{352, -44}, {369, -44}, {369, -42}, {384, -42}}, color = {0, 0, 127}));
  connect(omegaPu, add2.y) annotation(
    Line(points = {{414, 40}, {286, 40}}, color = {0, 0, 127}));
  connect(integrator1.y, phaseurGrid.UPhaseOffs) annotation(
    Line(points = {{408, -42}, {504, -42}, {504, -32}}, color = {0, 0, 127}));
  connect(Wref.y, add4.u2) annotation(
    Line(points = {{238, -8}, {252, -8}, {252, -50}, {328, -50}}, color = {0, 0, 127}));
  connect(Wref.y, add2.u2) annotation(
    Line(points = {{238, -8}, {244, -8}, {244, 34}, {262, 34}}, color = {0, 0, 127}));
  connect(phaseurGrid.PPu, add3.u2) annotation(
    Line(points = {{574, -4}, {634, -4}, {634, -104}, {50, -104}, {50, 0}, {72, 0}, {72, 50}, {94, 50}}, color = {0, 0, 127}));
  connect(phaseurGrid.terminal, aCPower) annotation(
    Line(points = {{572, -16}, {750, -16}}, color = {0, 0, 255}));
  connect(gain.u, gain3.y) annotation(
    Line(points = {{-8, -34}, {116.5, -34}, {116.5, -35}}, color = {0, 0, 127}));
  connect(gain3.u, integrator.y) annotation(
    Line(points = {{151, -35}, {174, -35}, {174, 50}}, color = {0, 0, 127}));
  annotation(
    Diagram(coordinateSystem(extent = {{-200, 100}, {780, -140}})),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002),
  Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}), graphics = {Text(origin = {175, -38}, extent = {{-45, 40}, {45, -40}}, textString = "OmegaPu"), Rectangle(extent = {{-100, 100}, {100, -100}}), Text(origin = {2, 8}, extent = {{-74, 50}, {74, -50}}, textString = "ACGrid")}));
end AcGrid;