#include "../../include/nuchic/settings.hh"

SettingsNode::SettingsNode() {
    parent = nullptr;
    children.resize(0);
}

SettingsNode::SettingsNode(SettingsPtr _parent) : parent(_parent) {
    children.resize(0);
}
