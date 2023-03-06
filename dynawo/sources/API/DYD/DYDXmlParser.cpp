//
// Copyright (c) 2023, RTE (http://www.rte-france.com)
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
 * @file DYDXmlParser.cpp
 * @brief XML parser : implementation file
 *
 */

#include "DYDXmlParser.h"
#include "DYNCommon.h"
#include "DYNXmlString.h"
#include "DYNXmlCharConversion.h"
#include "DYDModelicaModel.h"
#include "DYDModelTemplate.h"
#include "DYDBlackBoxModel.h"
#include "DYDModelTemplateExpansion.h"
#include "DYNMacrosMessage.h"
#include "DYDConnector.h"

#pragma GCC diagnostic ignored "-Wunused-parameter"

namespace dynamicdata {

void validationErrorHandler(void* /*userData*/, xmlErrorPtr error) {
    std::string message = "XML error: ";
    if (error != nullptr) {
        if (error->message != nullptr) {
            message += error->message;
        }
        if (error->file != nullptr) {
            message += " in file ";
            message += error->file;
        }
        message += " at line ";
        message += std::to_string(error->line);
    }
    throw DYNError(DYN::Error::API, XsdValidationError, message);
}

XmlParser::XmlParser(boost::shared_ptr<DynamicModelsCollection> dynamicModelsCollection, const std::string& filename) :
    dynamicModelsCollection_(dynamicModelsCollection),
    reader_(xmlNewTextReaderFilename(filename.c_str())),
    filename_(filename) {
    if (!reader_) {
        throw DYNError(DYN::Error::API, FileSystemItemDoesNotExist, filename);
    }
    xmlSetStructuredErrorFunc(NULL, &validationErrorHandler);
}

void
XmlParser::parseXML() {
    try {
        while (xmlTextReaderRead(reader_.get()) == 1) {
            if (xmlTextReaderNodeType(reader_.get()) == XML_READER_TYPE_ELEMENT) {
                const DYN::XmlString nodeName(xmlTextReaderName(reader_.get()));
                if (xmlStrEqual(nodeName.get(), DYN::S2XML("dyn:dynamicModelsArchitecture"))) {
                    parseDynamicModelsArchitecture();
                }
            }
        }
    } catch (const std::exception& err) {
        throw DYNError(DYN::Error::API, XmlFileParsingError, filename_, err.what());
    }
}

void
XmlParser::activateXSDValidation(const std::string& xsdSchemaFile) {
    xmlTextReaderSchemaValidate(reader_.get(), xsdSchemaFile.c_str());
}

void
XmlParser::parseDynamicModelsArchitecture() {
    while (xmlTextReaderRead(reader_.get()) == 1) {
        if (xmlTextReaderNodeType(reader_.get()) == XML_READER_TYPE_ELEMENT) {
            const DYN::XmlString nodeName(xmlTextReaderName(reader_.get()));
            if (xmlStrEqual(nodeName.get(), DYN::S2XML("dyn:modelicaModel"))) {
                parseModelicaModel();
            } else if (xmlStrEqual(nodeName.get(), DYN::S2XML("dyn:modelTemplate"))) {
                parseModelTemplate();
            } else if (xmlStrEqual(nodeName.get(), DYN::S2XML("dyn:blackBoxModel"))) {
                parseBlackBoxModel();
            } else if (xmlStrEqual(nodeName.get(), DYN::S2XML("dyn:modelTemplateExpansion"))) {
                parseModelTemplateExpansion();
            } else if (xmlStrEqual(nodeName.get(), DYN::S2XML("dyn:connect"))) {
                boost::shared_ptr<Connector> connector = parseConnect();
                dynamicModelsCollection_->addConnect(connector);
            } else if (xmlStrEqual(nodeName.get(), DYN::S2XML("dyn:macroConnect"))) {
                parseMacroConnect();
            } else if (xmlStrEqual(nodeName.get(), DYN::S2XML("dyn:macroConnector"))) {
                parseMacroConnector();
            } else if (xmlStrEqual(nodeName.get(), DYN::S2XML("dyn:macroStaticReference"))) {
                parseMacroStaticReference();
            } else {
                throw DYNError(DYN::Error::API, XmlUnknownNodeName, nodeName.get(), filename_);
            }
        } else if (xmlTextReaderNodeType(reader_.get()) == XML_READER_TYPE_END_ELEMENT) {
            const DYN::XmlString nodeName(xmlTextReaderName(reader_.get()));
            if (xmlStrEqual(nodeName.get(), DYN::S2XML("dyn:dynamicModelsArchitecture"))) {
                break;
            }
        }
    }
}

void
XmlParser::parseModelicaModel() {
    try {
        boost::shared_ptr<ModelicaModel> modelicaModel = nullptr;
        const DYN::XmlString modelicaModelId(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("id")));
        if (modelicaModelId != nullptr) {
            modelicaModel = boost::shared_ptr<ModelicaModel>(new ModelicaModel(DYN::XML2S(modelicaModelId.get())));
            const DYN::XmlString modelicaModelStaticId(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("staticId")));
            const DYN::XmlString modelicaModelUseAliasing(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("useAliasing")));
            const DYN::XmlString modelicaModelGenerateCalculatedVariables(
                                xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("generateCalculatedVariables")));
            if (modelicaModelStaticId != nullptr) {
                modelicaModel->setStaticId(DYN::XML2S(modelicaModelStaticId.get()));
            }
            if (modelicaModelUseAliasing != nullptr && modelicaModelGenerateCalculatedVariables != nullptr) {
                modelicaModel->setCompilationOptions(DYN::str2Bool(DYN::XML2S(modelicaModelUseAliasing.get())),
                                                        DYN::str2Bool(DYN::XML2S(modelicaModelGenerateCalculatedVariables.get())));
            }
        }

        while (xmlTextReaderRead(reader_.get()) == 1 && xmlTextReaderNodeType(reader_.get()) != XML_READER_TYPE_END_ELEMENT) {
            if (xmlTextReaderNodeType(reader_.get()) == XML_READER_TYPE_ELEMENT) {
                try {
                    const DYN::XmlString nodeName(xmlTextReaderName(reader_.get()));
                    if (xmlStrEqual(nodeName.get(), DYN::S2XML("dyn:unitDynamicModel"))) {
                        boost::shared_ptr<UnitDynamicModel> unitDynamicModel = parseUnitDynamicModel();
                        modelicaModel->addUnitDynamicModel(unitDynamicModel);
                    } else if (xmlStrEqual(nodeName.get(), DYN::S2XML("dyn:connect"))) {
                        boost::shared_ptr<Connector> connector = parseConnect();
                        modelicaModel->addConnect(connector);
                    } else if (xmlStrEqual(nodeName.get(), DYN::S2XML("dyn:initConnect"))) {
                        boost::shared_ptr<Connector> initConnector = parseConnect();
                        modelicaModel->addInitConnect(initConnector);
                    } else if (xmlStrEqual(nodeName.get(), DYN::S2XML("dyn::macroConnect"))) {
                        parseMacroConnect();  // pb
                    } else if (xmlStrEqual(nodeName.get(), DYN::S2XML("dyn:staticRef"))) {
                        parseStaticRef(modelicaModel);
                    } else if (xmlStrEqual(nodeName.get(), DYN::S2XML("dyn:macroStaticRef"))) {
                        parseMacroStaticRef(modelicaModel);
                    } else {
                        throw DYNError(DYN::Error::API, XmlUnknownNodeName, nodeName.get(), filename_);
                    }
                } catch (const std::exception& err) {
                    throw DYNError(DYN::Error::API, XmlFileParsingError, filename_, err.what());
                }
            }
        }
        dynamicModelsCollection_->addModel(modelicaModel);
    } catch (const std::exception& err) {
        throw DYNError(DYN::Error::API, XmlFileParsingError, filename_, err.what());
    }
}

void
XmlParser::parseModelTemplate() {
    // à remplir
    throw std::logic_error("À remplir");
}

void
XmlParser::parseBlackBoxModel() {
    try {
        boost::shared_ptr<BlackBoxModel> blackBoxModel = nullptr;
        const DYN::XmlString blackBoxModelId(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("id")));
        if (blackBoxModelId != nullptr) {
            blackBoxModel = boost::shared_ptr<BlackBoxModel>(new BlackBoxModel(DYN::XML2S(blackBoxModelId.get())));
            const DYN::XmlString blackBoxModelStaticId(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("staticId")));
            const DYN::XmlString blackBoxModelLib(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("lib")));
            const DYN::XmlString blackBoxModelParFile(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("parFile")));
            const DYN::XmlString blackBoxModelParId(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("parId")));
            if (blackBoxModelStaticId != nullptr) {
                blackBoxModel->setStaticId(DYN::XML2S(blackBoxModelStaticId.get()));
            }
            if (blackBoxModelLib != nullptr) {
                blackBoxModel->setLib(DYN::XML2S(blackBoxModelLib.get()));
            }
            if (blackBoxModelParFile != nullptr) {
                blackBoxModel->setParFile(DYN::XML2S(blackBoxModelParFile.get()));
            }
            if (blackBoxModelParId != nullptr) {
                blackBoxModel->setParId(DYN::XML2S(blackBoxModelParId.get()));
            }
        }
        parseBlackBoxModelOrModelTemplateExpansion(blackBoxModel);
        dynamicModelsCollection_->addModel(blackBoxModel);
    } catch (const std::exception& err) {
        throw DYNError(DYN::Error::API, XmlFileParsingError, filename_, err.what());
    }
}

void
XmlParser::parseModelTemplateExpansion() {
    try {
        boost::shared_ptr<ModelTemplateExpansion> modelTemplateExpansion = nullptr;
        const DYN::XmlString modelTemplateExpansionId(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("id")));
        if (modelTemplateExpansionId != nullptr) {
            modelTemplateExpansion = boost::shared_ptr<ModelTemplateExpansion>(
                                                    new ModelTemplateExpansion(DYN::XML2S(modelTemplateExpansionId.get())));
            const DYN::XmlString modelTemplateExpansionTemplateId(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("templateId")));
            const DYN::XmlString modelTemplateExpansionStaticId(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("staticId")));
            const DYN::XmlString modelTemplateExpansionParFile(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("parFile")));
            const DYN::XmlString modelTemplateExpansionParId(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("parId")));

            if (modelTemplateExpansionTemplateId != nullptr) {
                modelTemplateExpansion->setTemplateId(DYN::XML2S(modelTemplateExpansionTemplateId.get()));
            }
            if (modelTemplateExpansionStaticId != nullptr) {
                modelTemplateExpansion->setStaticId(DYN::XML2S(modelTemplateExpansionStaticId.get()));
            }
            if (modelTemplateExpansionParFile != nullptr) {
                modelTemplateExpansion->setParFile(DYN::XML2S(modelTemplateExpansionParFile.get()));
            }
            if (modelTemplateExpansionParId != nullptr) {
                modelTemplateExpansion->setParId(DYN::XML2S(modelTemplateExpansionParId.get()));
            }
        }
        parseBlackBoxModelOrModelTemplateExpansion(modelTemplateExpansion);
        dynamicModelsCollection_->addModel(modelTemplateExpansion);
    } catch (const std::exception& err) {
        throw DYNError(DYN::Error::API, XmlFileParsingError, filename_, err.what());
    }
}

void
XmlParser::parseMacroConnect() {
    // à remplir
    throw std::logic_error("À remplir");
}

void
XmlParser::parseMacroConnector() {
    // à remplir
    throw std::logic_error("À remplir");
}

void
XmlParser::parseMacroStaticReference() {
    // à remplir
    throw std::logic_error("À remplir");
}

boost::shared_ptr<Connector>
XmlParser::parseConnect() const {
    try {
        boost::shared_ptr<Connector> connector = nullptr;
        const DYN::XmlString connectorId1(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("id1")));
        const DYN::XmlString connectorVar1(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("var1")));
        const DYN::XmlString connectorId2(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("id2")));
        const DYN::XmlString connectorVar2(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("var2")));
        if (connectorId1 != nullptr && connectorVar1 != nullptr && connectorId2 != nullptr && connectorVar2 != nullptr) {
            connector = boost::shared_ptr<Connector>(new Connector(DYN::XML2S(connectorId1.get()),
                                                                    DYN::XML2S(connectorVar1.get()),
                                                                    DYN::XML2S(connectorId2.get()),
                                                                    DYN::XML2S(connectorVar2.get())));
        } else {
            throw DYNError(DYN::Error::API, XmlMissingAttributeError, filename_, "dyn:connect/dyn:initConnect");
        }
        return connector;
    } catch (const std::exception& err) {
        throw DYNError(DYN::Error::API, XmlFileParsingError, filename_, err.what());
    }
}

void
XmlParser::parseMacroConnection() const {
    // à remplir
    throw std::logic_error("À remplir");
}

boost::shared_ptr<UnitDynamicModel>
XmlParser::parseUnitDynamicModel() const {
    try {
        boost::shared_ptr<UnitDynamicModel> unitDynamicModel = nullptr;

        const DYN::XmlString unitDynamicModelId(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("id")));
        const DYN::XmlString unitDynamicModelDynamicModelName(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("name")));
        const DYN::XmlString unitDynamicModelDynamicFileName(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("moFile")));
        const DYN::XmlString unitDynamicModelInitModelName(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("initName")));
        const DYN::XmlString unitDynamicModelInitFileName(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("initFile")));
        const DYN::XmlString unitDynamicModelParFile(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("parFile")));
        const DYN::XmlString unitDynamicModelParId(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("parId")));

        if (unitDynamicModelId != nullptr && unitDynamicModelDynamicModelName != nullptr) {
            unitDynamicModel = boost::shared_ptr<UnitDynamicModel>(new UnitDynamicModel(DYN::XML2S(unitDynamicModelId.get()),
                                                                            DYN::XML2S(unitDynamicModelDynamicModelName.get())));
            if (unitDynamicModelDynamicFileName != nullptr) {
                unitDynamicModel->setDynamicFileName(DYN::XML2S(unitDynamicModelDynamicFileName.get()));
            }
            if (unitDynamicModelInitModelName != nullptr) {
                unitDynamicModel->setInitModelName(DYN::XML2S(unitDynamicModelInitModelName.get()));
            }
            if (unitDynamicModelInitFileName != nullptr) {
                unitDynamicModel->setInitFileName(DYN::XML2S(unitDynamicModelInitFileName.get()));
            }
            if (unitDynamicModelParFile != nullptr) {
                unitDynamicModel->setParFile(DYN::XML2S(unitDynamicModelParFile.get()));
            }
            if (unitDynamicModelParId != nullptr) {
                unitDynamicModel->setParId(DYN::XML2S(unitDynamicModelParId.get()));
            }
        } else {
            throw DYNError(DYN::Error::API, XmlMissingAttributeError, filename_, "dyn:unitDynamicModel");
        }

        return unitDynamicModel;
    } catch (const std::exception& err) {
        throw DYNError(DYN::Error::API, XmlFileParsingError, filename_, err.what());
    }
}

void
XmlParser::parseStaticRef(boost::shared_ptr<Model> model) const {
    try {
        model->addStaticRef(DYN::XML2S(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("var"))),
                            DYN::XML2S(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("staticVar"))));
    } catch (const std::exception& err) {
        throw DYNError(DYN::Error::API, XmlFileParsingError, filename_, err.what());
    }
}

void
XmlParser::parseMacroStaticRef(boost::shared_ptr<Model> model) const {
    try {
        boost::shared_ptr<MacroStaticRef> macroStaticRef = boost::shared_ptr<MacroStaticRef>(
                            new MacroStaticRef(DYN::XML2S(xmlTextReaderGetAttribute(reader_.get(), DYN::S2XML("id")))));
        model->addMacroStaticRef(macroStaticRef);
    } catch (const std::exception& err) {
        throw DYNError(DYN::Error::API, XmlFileParsingError, filename_, err.what());
    }
}

void
XmlParser::parseBlackBoxModelOrModelTemplateExpansion(boost::shared_ptr<Model> model) const {
    while (xmlTextReaderRead(reader_.get()) == 1 && xmlTextReaderNodeType(reader_.get()) != XML_READER_TYPE_END_ELEMENT) {
        if (xmlTextReaderNodeType(reader_.get()) == XML_READER_TYPE_ELEMENT) {
            try {
                const DYN::XmlString nodeName(xmlTextReaderName(reader_.get()));
                if (xmlStrEqual(nodeName.get(), DYN::S2XML("dyn:staticRef"))) {
                    parseStaticRef(model);
                } else if (xmlStrEqual(nodeName.get(), DYN::S2XML("dyn:macroStaticRef"))) {
                    parseMacroStaticRef(model);
                } else {
                    throw DYNError(DYN::Error::API, XmlUnknownNodeName, nodeName.get(), filename_);
                }
            } catch (const std::exception& err) {
                throw DYNError(DYN::Error::API, XmlFileParsingError, filename_, err.what());
            }
        }
    }
}

}  // namespace dynamicdata
