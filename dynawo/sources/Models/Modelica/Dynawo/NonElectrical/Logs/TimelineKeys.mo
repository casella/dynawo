within Dynawo.NonElectrical.Logs;

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

encapsulated package TimelineKeys

  final constant Integer ActivatePMAX = 0;
  final constant Integer ActivatePMIN = 1;
  final constant Integer BusAboveVoltage = 2;
  final constant Integer BusUnderVoltage = 3;
  final constant Integer ComponentConnected = 4;
  final constant Integer ComponentDisconnected = 5;
  final constant Integer Converter1Connected = 6;
  final constant Integer Converter1SwitchOff = 7;
  final constant Integer Converter2Connected = 8;
  final constant Integer Converter2SwitchOff = 9;
  final constant Integer CriteriaNotChecked = 10;
  final constant Integer CurrentLimitAutomatonActing = 11;
  final constant Integer CurrentLimitAutomatonArming = 12;
  final constant Integer CurrentLimitAutomatonDisarming = 13;
  final constant Integer DCLineClosed = 14;
  final constant Integer DCLineOpen = 15;
  final constant Integer DanglingLineConnected = 16;
  final constant Integer DanglingLineDisconnected = 17;
  final constant Integer DeactivatePMAX = 18;
  final constant Integer DeactivatePMIN = 19;
  final constant Integer DistanceTrippedZone1 = 20;
  final constant Integer DistanceTrippedZone2 = 21;
  final constant Integer DistanceTrippedZone3 = 22;
  final constant Integer DistanceTrippedZone4 = 23;
  final constant Integer GeneratorBackRegulation = 24;
  final constant Integer GeneratorConnected = 25;
  final constant Integer GeneratorDisconnected = 26;
  final constant Integer GeneratorMaxQ = 27;
  final constant Integer GeneratorMinQ = 28;
  final constant Integer GeneratorPVBackRegulation = 29;
  final constant Integer GeneratorPVMaxQ = 30;
  final constant Integer GeneratorPVMinQ = 31;
  final constant Integer GeneratorTargetP = 32;
  final constant Integer GeneratorTargetQ = 33;
  final constant Integer HVDC1BackRegulation = 34;
  final constant Integer HVDC1MaxQ = 35;
  final constant Integer HVDC1MinQ = 36;
  final constant Integer HVDC2BackRegulation = 37;
  final constant Integer HVDC2MaxQ = 38;
  final constant Integer HVDC2MinQ = 39;
  final constant Integer HVDCDeactivateMaxP = 40;
  final constant Integer HVDCMaxP = 41;
  final constant Integer IdealSwitchSwitchOff = 42;
  final constant Integer IdealSwitchSwitchOn = 43;
  final constant Integer LVRTArming = 44;
  final constant Integer LVRTDisarming = 45;
  final constant Integer LVRTTripped = 46;
  final constant Integer LineCloseSide1 = 47;
  final constant Integer LineCloseSide2 = 48;
  final constant Integer LineClosed = 49;
  final constant Integer LineOpen = 50;
  final constant Integer LineOpenSide1 = 51;
  final constant Integer LineOpenSide2 = 52;
  final constant Integer LoadConnected = 53;
  final constant Integer LoadDisconnected = 54;
  final constant Integer LoadModificationEnded = 55;
  final constant Integer LoadModificationStarted = 56;
  final constant Integer LoadSheddingStarted = 57;
  final constant Integer LossOfSynchronismArming = 58;
  final constant Integer LossOfSynchronismTripped = 59;
  final constant Integer NodeFaultBegin = 60;
  final constant Integer NodeFaultEnd = 61;
  final constant Integer NodeOff = 62;
  final constant Integer NodeOn = 63;
  final constant Integer OverloadDown = 64;
  final constant Integer OverloadOpen = 65;
  final constant Integer OverloadUp = 66;
  final constant Integer OverspeedArming = 67;
  final constant Integer OverspeedDisarming = 68;
  final constant Integer OverspeedTripped = 69;
  final constant Integer PhaseShifterAboveMax = 70;
  final constant Integer PhaseShifterBelowMax = 71;
  final constant Integer PhaseShifterBelowMin = 72;
  final constant Integer PhaseShifterBelowStop = 73;
  final constant Integer PhaseShifterBlockingIActing = 74;
  final constant Integer PhaseShifterBlockingIArming = 75;
  final constant Integer PhaseShifterBlockingIDisarming = 76;
  final constant Integer PhaseShifterSwitchOff = 77;
  final constant Integer PhaseShifterSwitchOn = 78;
  final constant Integer PhaseShifterWithinInterval = 79;
  final constant Integer ProtectionConnected = 80;
  final constant Integer ProtectionDisconnected = 81;
  final constant Integer RPCLLimitationUsRefMax = 82;
  final constant Integer RPCLLimitationUsRefMin = 83;
  final constant Integer RPCLStandard = 84;
  final constant Integer SVarCBackRegulation = 85;
  final constant Integer SVarCConnected = 86;
  final constant Integer SVarCDisconnected = 87;
  final constant Integer SVarCMaxB = 88;
  final constant Integer SVarCMinB = 89;
  final constant Integer SVarCOff = 90;
  final constant Integer SVarCRunning = 91;
  final constant Integer SVarCStandby = 92;
  final constant Integer SVarCUmaxreached = 93;
  final constant Integer SVarCUminreached = 94;
  final constant Integer ShuntConnected = 95;
  final constant Integer ShuntDisconnected = 96;
  final constant Integer SignalReceived = 97;
  final constant Integer SourceAbovePower = 98;
  final constant Integer SourcePowerAboveMax = 99;
  final constant Integer SourcePowerBelowMin = 100;
  final constant Integer SourcePowerTakenIntoAccount = 101;
  final constant Integer SourceUnderPower = 102;
  final constant Integer SwitchClosed = 103;
  final constant Integer SwitchOpened = 104;
  final constant Integer TapChangerAboveMax = 105;
  final constant Integer TapChangerBelowMin = 106;
  final constant Integer TapChangerSwitchOff = 107;
  final constant Integer TapChangerSwitchOn = 108;
  final constant Integer TapChangersArming = 109;
  final constant Integer TapChangersBlocked = 110;
  final constant Integer TapChangersBlockedD = 111;
  final constant Integer TapChangersBlockedT = 112;
  final constant Integer TapChangersUnarming = 113;
  final constant Integer TapChangersUnblocked = 114;
  final constant Integer TapDown = 115;
  final constant Integer TapUp = 116;
  final constant Integer TerminateInModel = 117;
  final constant Integer TransformerSwitchOff = 118;
  final constant Integer TransformerSwitchOn = 119;
  final constant Integer TwoWTFOCloseSide1 = 120;
  final constant Integer TwoWTFOCloseSide2 = 121;
  final constant Integer TwoWTFOClosed = 122;
  final constant Integer TwoWTFOOpen = 123;
  final constant Integer TwoWTFOOpenSide1 = 124;
  final constant Integer TwoWTFOOpenSide2 = 125;
  final constant Integer UFLS10Activated = 126;
  final constant Integer UFLS10Arming = 127;
  final constant Integer UFLS1Activated = 128;
  final constant Integer UFLS1Arming = 129;
  final constant Integer UFLS2Activated = 130;
  final constant Integer UFLS2Arming = 131;
  final constant Integer UFLS3Activated = 132;
  final constant Integer UFLS3Arming = 133;
  final constant Integer UFLS4Activated = 134;
  final constant Integer UFLS4Arming = 135;
  final constant Integer UFLS5Activated = 136;
  final constant Integer UFLS5Arming = 137;
  final constant Integer UFLS6Activated = 138;
  final constant Integer UFLS6Arming = 139;
  final constant Integer UFLS7Activated = 140;
  final constant Integer UFLS7Arming = 141;
  final constant Integer UFLS8Activated = 142;
  final constant Integer UFLS8Arming = 143;
  final constant Integer UFLS9Activated = 144;
  final constant Integer UFLS9Arming = 145;
  final constant Integer UVAArming = 146;
  final constant Integer UVADisarming = 147;
  final constant Integer UVATripped = 148;
  final constant Integer UnderspeedArming = 149;
  final constant Integer UnderspeedDisarming = 150;
  final constant Integer UnderspeedTripped = 151;
  final constant Integer VRBackToRegulation = 152;
  final constant Integer VRFrozen = 153;
  final constant Integer VRLimitationEfdMax = 154;
  final constant Integer VRLimitationEfdMin = 155;
  final constant Integer VRLimitationUsRefMax = 156;
  final constant Integer VRLimitationUsRefMin = 157;
  final constant Integer VRUnfrozen = 158;
  final constant Integer VoltageSetPointChangeEnded = 159;
  final constant Integer VoltageSetPointChangeStarted = 160;
  final constant Integer Zone1Arming = 161;
  final constant Integer Zone1Disarming = 162;
  final constant Integer Zone2Arming = 163;
  final constant Integer Zone2Disarming = 164;
  final constant Integer Zone3Arming = 165;
  final constant Integer Zone3Disarming = 166;
  final constant Integer Zone4Arming = 167;
  final constant Integer Zone4Disarming = 168;

  annotation(preferredView = "text");
end TimelineKeys;
