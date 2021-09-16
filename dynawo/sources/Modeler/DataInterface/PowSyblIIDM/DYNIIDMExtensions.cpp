//
// Copyright (c) 2021, RTE (http://www.rte-france.com)
// See AUTHORS.txt
// All rights reserved.
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
//
// This file is part of Dynawo, an hybrid C++/Modelica open source time domain
// simulation tool for power systems.
//

/**
 * @file DYNIIDMExtensions.cpp
 * @brief File for external IIDM extensions management : implementation
 */

#include "DYNIIDMExtensions.hpp"

#include "DYNExecUtils.h"
#include "DYNFileSystemUtils.h"
#include "DYNMacrosMessage.h"

namespace DYN {

std::unordered_map<IIDMExtensions::LibraryPath, std::shared_ptr<boost::dll::shared_library> > IIDMExtensions::libraries_;
std::unordered_set<IIDMExtensions::LibraryPath> IIDMExtensions::librariesLoadingIssues_;
std::mutex IIDMExtensions::librariesMutex_;

boost::filesystem::path
IIDMExtensions::findLibraryPath() {
  auto libPathVar = getenv("DYNAWO_IIDM_EXTENSION");
  if (!libPathVar) {
    return boost::filesystem::path();
  }
  auto libPath = std::string(libPathVar);
  if (!exists(libPath)) {
    return boost::filesystem::path();
  }
  return boost::filesystem::path(libPath);
}

}  // namespace DYN