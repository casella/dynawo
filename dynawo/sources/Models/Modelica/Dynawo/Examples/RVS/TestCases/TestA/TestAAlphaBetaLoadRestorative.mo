within Dynawo.Examples.RVS.TestCases.TestA;

/*
* Copyright (c) 2023, RTE (http://www.rte-france.com)
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

model TestAAlphaBetaLoadRestorative "RVS test system simulation case : reactive load connection, restorative loads"
  import Modelica.SIunits.Conversions.from_deg;

  extends Icons.Example;
  extends Dynawo.Examples.RVS.Grid.FullDynamicInfiniteBus(
    P0Pu_gen_10101_ABEL_G1 = -0.09999999999492322,
    P0Pu_gen_10102_ADAMS_G1 = -0.09999999999494089,
    P0Pu_gen_10107_ALDER_G1 = -0.7999999999726052,
    P0Pu_gen_10113_ARNE_G1 = -1.6250000001927132,
    P0Pu_gen_10115_ARTHUR_G1 = -0.11999999999728818,
    P0Pu_gen_10116_ASSER_G1 = -1.550000148856329,
    P0Pu_gen_10118_ASTOR_G1 = -4.000000355338675,
    P0Pu_gen_10122_AUBREY_G1 = -0.5000004099960184,
    P0Pu_gen_10123_AUSTEN_G1 = -1.550000063473528,
    P0Pu_gen_20101_ABEL_G2 = -0.09999999999492801,
    P0Pu_gen_20102_ADAMS_G2 = -0.09999999999493148,
    P0Pu_gen_20107_ALDER_G2 = -0.7999999999725957,
    P0Pu_gen_20113_ARNE_G2 = -1.6249999999306455,
    P0Pu_gen_20115_ARTHUR_G2 = -0.11999999999728826,
    P0Pu_gen_20122_AUBREY_G2 = -0.5000004099960331,
    P0Pu_gen_20123_AUSTEN_G2 = -1.5500000634748163,
    P0Pu_gen_30101_ABEL_G3 = -0.7599999999788776,
    P0Pu_gen_30102_ADAMS_G3 = -0.7599999999782525,
    P0Pu_gen_30107_ALDER_G3 = -0.7999999999726064,
    P0Pu_gen_30113_ARNE_G3 = -1.6249999999307339,
    P0Pu_gen_30115_ARTHUR_G3 = -0.11999999999728914,
    P0Pu_gen_30122_AUBREY_G3 = -0.5000004099960347,
    P0Pu_gen_30123_AUSTEN_G3 = -3.5000001656289816,
    P0Pu_gen_40101_ABEL_G4 = -0.7599999999788802,
    P0Pu_gen_40102_ADAMS_G4 = -0.7599999999782637,
    P0Pu_gen_40115_ARTHUR_G4 = -0.1199999999972882,
    P0Pu_gen_40122_AUBREY_G4 = -0.5000004099960342,
    P0Pu_gen_50115_ARTHUR_G5 = -0.1199999999972885,
    P0Pu_gen_50122_AUBREY_G5 = -0.5000004099960343,
    P0Pu_gen_60115_ARTHUR_G6 = -1.550000141046319,
    P0Pu_gen_60122_AUBREY_G6 = -0.5000004099960347,
    Q0Pu_gen_10101_ABEL_G1 = -0.0747006931432781,
    Q0Pu_gen_10102_ADAMS_G1 = -0.10294511872209362,
    Q0Pu_gen_10107_ALDER_G1 = -0.5106708523802469,
    Q0Pu_gen_10113_ARNE_G1 = -0.0378376197136411,
    Q0Pu_gen_10115_ARTHUR_G1 = -0.06383030028772678,
    Q0Pu_gen_10116_ASSER_G1 = -0.5361276776730247,
    Q0Pu_gen_10118_ASTOR_G1 = -2.0301805325991826,
    Q0Pu_gen_10122_AUBREY_G1 = -0.036463694362460664,
    Q0Pu_gen_10123_AUSTEN_G1 = -0.7087170883544405,
    Q0Pu_gen_20101_ABEL_G2 = -0.07470069314325702,
    Q0Pu_gen_20102_ADAMS_G2 = -0.10294511872206842,
    Q0Pu_gen_20107_ALDER_G2 = -0.5107086749439578,
    Q0Pu_gen_20113_ARNE_G2 = -1.2008084007899797,
    Q0Pu_gen_20115_ARTHUR_G2 = -0.06383030028773819,
    Q0Pu_gen_20122_AUBREY_G2 = -0.03646369436252084,
    Q0Pu_gen_20123_AUSTEN_G2 = -0.7087170883555296,
    Q0Pu_gen_30101_ABEL_G3 = -0.22178314963648094,
    Q0Pu_gen_30102_ADAMS_G3 = -0.18651193763838703,
    Q0Pu_gen_30107_ALDER_G3 = -0.5106708523801484,
    Q0Pu_gen_30113_ARNE_G3 = -1.2008829225627582,
    Q0Pu_gen_30115_ARTHUR_G3 = -0.0638303002877329,
    Q0Pu_gen_30122_AUBREY_G3 = -0.03646369436251523,
    Q0Pu_gen_30123_AUSTEN_G3 = -1.3305380937734208,
    Q0Pu_gen_40101_ABEL_G4 = -0.22175966597865165,
    Q0Pu_gen_40102_ADAMS_G4 = -0.18651193763830845,
    Q0Pu_gen_40115_ARTHUR_G4 = -0.06383987947986561,
    Q0Pu_gen_40122_AUBREY_G4 = -0.03646369436253674,
    Q0Pu_gen_50115_ARTHUR_G5 = -0.06383030028772911,
    Q0Pu_gen_50122_AUBREY_G5 = -0.0364636943624852,
    Q0Pu_gen_60115_ARTHUR_G6 = -0.8602969818617706,
    Q0Pu_gen_60122_AUBREY_G6 = -0.03646369436248653,
    U0Pu_gen_10101_ABEL_G1 = 1.0296020605827343,
    U0Pu_gen_10102_ADAMS_G1 = 1.047286860269937,
    U0Pu_gen_10107_ALDER_G1 = 1.0392093432573648,
    U0Pu_gen_10113_ARNE_G1 = 0.9699722451585717,
    U0Pu_gen_10115_ARTHUR_G1 = 1.0237450113009863,
    U0Pu_gen_10116_ASSER_G1 = 1.005707605152536,
    U0Pu_gen_10118_ASTOR_G1 = 1.049501524085332,
    U0Pu_gen_10122_AUBREY_G1 = 1.0025914159367653,
    U0Pu_gen_10123_AUSTEN_G1 = 1.0505080247133123,
    U0Pu_gen_20101_ABEL_G2 = 1.029602060582723,
    U0Pu_gen_20102_ADAMS_G2 = 1.047286860269925,
    U0Pu_gen_20107_ALDER_G2 = 1.0392137465540379,
    U0Pu_gen_20113_ARNE_G2 = 1.0427229897892534,
    U0Pu_gen_20115_ARTHUR_G2 = 1.023745011300991,
    U0Pu_gen_20122_AUBREY_G2 = 1.002591415936771,
    U0Pu_gen_20123_AUSTEN_G2 = 1.0505080247133731,
    U0Pu_gen_30101_ABEL_G3 = 1.0161596964727235,
    U0Pu_gen_30102_ADAMS_G3 = 1.0119115875875173,
    U0Pu_gen_30107_ALDER_G3 = 1.0392093432573528,
    U0Pu_gen_30113_ARNE_G3 = 1.042727343363208,
    U0Pu_gen_30115_ARTHUR_G3 = 1.023745011300996,
    U0Pu_gen_30122_AUBREY_G3 = 1.0025914159367775,
    U0Pu_gen_30123_AUSTEN_G3 = 1.041325051732186,
    U0Pu_gen_40101_ABEL_G4 = 1.0161558787590688,
    U0Pu_gen_40102_ADAMS_G4 = 1.0119115875874976,
    U0Pu_gen_40115_ARTHUR_G4 = 1.023754566693637,
    U0Pu_gen_40122_AUBREY_G4 = 1.002591415936775,
    U0Pu_gen_50115_ARTHUR_G5 = 1.0237450113009892,
    U0Pu_gen_50122_AUBREY_G5 = 1.0025914159367704,
    U0Pu_gen_60115_ARTHUR_G6 = 1.0261638226400809,
    U0Pu_gen_60122_AUBREY_G6 = 1.00259141593677,
    UPhase0_gen_10101_ABEL_G1 = -0.26528375249349523,
    UPhase0_gen_10102_ADAMS_G1 = -0.26873055958340697,
    UPhase0_gen_10107_ALDER_G1 = -0.2909981796235494,
    UPhase0_gen_10113_ARNE_G1 = 0.00410841273654951,
    UPhase0_gen_10115_ARTHUR_G1 = 0.14591107650805835,
    UPhase0_gen_10116_ASSER_G1 = 0.13728705789917342,
    UPhase0_gen_10118_ASTOR_G1 = 0.219179365403202,
    UPhase0_gen_10122_AUBREY_G1 = 0.35407238234960037,
    UPhase0_gen_10123_AUSTEN_G1 = 0.15217939470869046,
    UPhase0_gen_20101_ABEL_G2 = -0.265283752493495,
    UPhase0_gen_20102_ADAMS_G2 = -0.2687305595834085,
    UPhase0_gen_20107_ALDER_G2 = -0.29099869448046284,
    UPhase0_gen_20113_ARNE_G2 = -0.005206386553900649,
    UPhase0_gen_20115_ARTHUR_G2 = 0.14591107650805607,
    UPhase0_gen_20122_AUBREY_G2 = 0.3540723823496088,
    UPhase0_gen_20123_AUSTEN_G2 = 0.15217939470878503,
    UPhase0_gen_30101_ABEL_G3 = -0.1984468778479227,
    UPhase0_gen_30102_ADAMS_G3 = -0.1999727128428949,
    UPhase0_gen_30107_ALDER_G3 = -0.29099817962354707,
    UPhase0_gen_30113_ARNE_G3 = -0.005206911147618021,
    UPhase0_gen_30115_ARTHUR_G3 = 0.14591107650805565,
    UPhase0_gen_30122_AUBREY_G3 = 0.35407238234960403,
    UPhase0_gen_30123_AUSTEN_G3 = 0.1531333025151461,
    UPhase0_gen_40101_ABEL_G4 = -0.19844631612764815,
    UPhase0_gen_40102_ADAMS_G4 = -0.19997271284288592,
    UPhase0_gen_40115_ARTHUR_G4 = 0.1459096513171526,
    UPhase0_gen_40122_AUBREY_G4 = 0.3540723823496089,
    UPhase0_gen_50115_ARTHUR_G5 = 0.14591107650805663,
    UPhase0_gen_50122_AUBREY_G5 = 0.35407238234960603,
    UPhase0_gen_60115_ARTHUR_G6 = 0.14471146857702907,
    UPhase0_gen_60122_AUBREY_G6 = 0.3540723823496069,
    Q0Pu_line_reactor_106 = 0.75,
    Q0Pu_line_reactor_110 = 0.75,
    U0Pu_line_reactor_106 = 1,
    U0Pu_line_reactor_110 = 1,
    UPhase0_line_reactor_106 = 0,
    UPhase0_line_reactor_110 = 0,
    P0Pu_load_1101_ABEL = 1.1880000305126062,
    P0Pu_load_1102_ADAMS = 1.0669999694032049,
    P0Pu_load_1103_ADLER = 1.9799999999965128,
    P0Pu_load_1104_AGRICOLA = 0.8140000152415523,
    P0Pu_load_1105_AIKEN = 0.7809999847285005,
    P0Pu_load_1106_ALBER = 1.4960000610362467,
    P0Pu_load_1107_ALDER = 1.3749999999783746,
    P0Pu_load_1108_ALGER = 1.8810000610320796,
    P0Pu_load_1109_ALI = 1.9249999999952647,
    P0Pu_load_1110_ALLEN = 2.144999999966342,
    P0Pu_load_1113_ARNE = 2.9149999999494955,
    P0Pu_load_1114_ARNOLD = 2.1339999389236377,
    P0Pu_load_1115_ARTHUR = 3.4870001220067244,
    P0Pu_load_1116_ASSER = 1.0999999999349863,
    P0Pu_load_1118_ASTOR = 3.662999877912528,
    P0Pu_load_1119_ATTAR = 1.9910000610193985,
    P0Pu_load_1120_ATTILA = 1.408000030459867,
    Q0Pu_load_1101_ABEL = 0.24200000762794033,
    Q0Pu_load_1102_ADAMS = 0.21999999997414654,
    Q0Pu_load_1103_ADLER = 0.407000007625965,
    Q0Pu_load_1104_AGRICOLA = 0.16499999998875836,
    Q0Pu_load_1105_AIKEN = 0.15399999617863952,
    Q0Pu_load_1106_ALBER = 0.3079999923711575,
    Q0Pu_load_1107_ALDER = 0.27499999999505487,
    Q0Pu_load_1108_ALGER = 0.3849999999997419,
    Q0Pu_load_1109_ALI = 0.3959999847378253,
    Q0Pu_load_1110_ALLEN = 0.4399999999884504,
    Q0Pu_load_1113_ARNE = 0.5940000152398809,
    Q0Pu_load_1114_ARNOLD = 0.42900001524297804,
    Q0Pu_load_1115_ARTHUR = 0.7040000152345325,
    Q0Pu_load_1116_ASSER = 0.21999999997646183,
    Q0Pu_load_1118_ASTOR = 0.7480000305135921,
    Q0Pu_load_1119_ATTAR = 0.40700000762424227,
    Q0Pu_load_1120_ATTILA = 0.286000003796252,
    U0Pu_load_1101_ABEL = 1.0429027751573068,
    U0Pu_load_1102_ADAMS = 1.041054780995125,
    U0Pu_load_1103_ADLER = 1.0477531875150454,
    U0Pu_load_1104_AGRICOLA = 1.0487827672892787,
    U0Pu_load_1105_AIKEN = 1.0473536193034214,
    U0Pu_load_1106_ALBER = 1.0354397985802706,
    U0Pu_load_1107_ALDER = 1.0422919000443116,
    U0Pu_load_1108_ALGER = 1.0434558747397327,
    U0Pu_load_1109_ALI = 1.048299224515768,
    U0Pu_load_1110_ALLEN = 1.0418891319680037,
    U0Pu_load_1113_ARNE = 1.047455268657455,
    U0Pu_load_1114_ARNOLD = 1.0414539630157793,
    U0Pu_load_1115_ARTHUR = 1.0429028244199405,
    U0Pu_load_1116_ASSER = 1.040889111806529,
    U0Pu_load_1118_ASTOR = 1.0434652629379533,
    U0Pu_load_1119_ATTAR = 1.0465844529211161,
    U0Pu_load_1120_ATTILA = 1.047790160113538,
    UPhase0_load_1101_ABEL = -0.43182718990200003,
    UPhase0_load_1102_ADAMS = -0.4236372139990933,
    UPhase0_load_1103_ADLER = -0.3893664931386952,
    UPhase0_load_1104_AGRICOLA = -0.46474061045111326,
    UPhase0_load_1105_AIKEN = -0.47039879672533297,
    UPhase0_load_1106_ALBER = -0.5170487955808892,
    UPhase0_load_1107_ALDER = -0.4821319549172829,
    UPhase0_load_1108_ALGER = -0.5018347629265315,
    UPhase0_load_1109_ALI = -0.4061940044819785,
    UPhase0_load_1110_ALLEN = -0.461408328490733,
    UPhase0_load_1113_ARNE = -0.21772124333877546,
    UPhase0_load_1114_ARNOLD = -0.28397712045663737,
    UPhase0_load_1115_ARTHUR = -0.09965058452637973,
    UPhase0_load_1116_ASSER = -0.09205886209200101,
    UPhase0_load_1118_ASTOR = -0.010656304687427738,
    UPhase0_load_1119_ATTAR = -0.11868130187462118,
    UPhase0_load_1120_ATTILA = -0.08659618500103589,
    P0Pu_sVarC_10106_ALBER_SVC = -4.0495384823202585e-14,
    P0Pu_sVarC_10114_ARNOLD_SVC = -8.97615315409439e-14,
    Q0Pu_sVarC_10106_ALBER_SVC = -0.6164618807059862,
    Q0Pu_sVarC_10114_ARNOLD_SVC = -1.254762105260428,
    U0Pu_sVarC_10106_ALBER_SVC = 1.049999952316467,
    U0Pu_sVarC_10114_ARNOLD_SVC = 1.049999952316385,
    UPhase0_sVarC_10106_ALBER_SVC = -0.4170880483313303,
    UPhase0_sVarC_10114_ARNOLD_SVC = -0.17175578394927937);

  parameter Real event_time_loadstep = 1;
  parameter Real event_QPu_load_101 = 2.5;
  parameter Real event_time_line_sw1 = 100;
  parameter Real event_time_line_sw2 = 100;
  parameter Real event_time_reactor110_sw1 = 100;
  parameter Real event_time_reactor110_sw2 = 100;
  parameter Real event_time_reactor106_sw1 = 100;
  parameter Real event_time_reactor106_sw2 = 100;

  Dynawo.Electrical.Loads.LoadAlphaBeta load_101(alpha = 1, beta = 2, i0Pu = Complex(0, 0), s0Pu = Complex(0, 0), u0Pu = ComplexMath.fromPolar(1.034, from_deg(-18.6))) annotation(
    Placement(visible = true, transformation(origin = {-310, -206}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));

equation
  line_106_110.switchOffSignal1.value = if time < event_time_line_sw1 then false else true;
  line_106_110.switchOffSignal2.value = if time < event_time_line_sw2 then false else true;
  line_reactor_110.switchOffSignal1.value = if time < event_time_reactor110_sw1 then false else true;
  line_reactor_110.switchOffSignal2.value = if time < event_time_reactor110_sw2 then false else true;
  line_reactor_106.switchOffSignal1.value = if time < event_time_reactor106_sw1 then false else true;
  line_reactor_106.switchOffSignal2.value = if time < event_time_reactor106_sw2 then false else true;
  load_101.deltaP = 0;
  load_101.deltaQ = 0;
  load_101.PRefPu = 0;
  load_101.QRefPu = if time < event_time_loadstep then 0 else event_QPu_load_101;
  load_101.switchOffSignal1.value = false;
  load_101.switchOffSignal2.value = false;

  connect(load_101.terminal, bus_101_ABEL.terminal) annotation(
    Line(points = {{-310, -206}, {-168, -206}}, color = {0, 0, 255}));

  annotation(preferredView = "diagram",
    experiment(StartTime = 0, StopTime = 80, Tolerance = 0.005, Interval = 0.01),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian --daeMode",
    __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "ida"),
    Diagram(coordinateSystem(extent = {{-340, -340}, {340, 340}})));
end TestAAlphaBetaLoadRestorative;
