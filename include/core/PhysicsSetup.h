#pragma once

#include <core/LQCDataTypes.h>
#include <IOdata/FileReader.h>
#include <IOdata/FastParser.h>

#include <filesystem>
#include <map>
#include <string>
#include <vector>
#include <stdexcept>
#include <print>

namespace lqcd {

/**
 * @brief 自动定位 docs/physics_setup.json 的绝对或相对路径
 */
[[nodiscard]] inline std::filesystem::path resolve_physics_setup_path(
    const std::filesystem::path& user_hint = "docs/physics_setup.json") {
    if (std::filesystem::exists(user_hint)) {
        return user_hint;
    }
    const std::filesystem::path candidates[] = {
        "docs/physics_setup.json",
        "../docs/physics_setup.json",
        "../../docs/physics_setup.json",
        "/Users/junxiongnie/code/build/ParseLQCData/docs/physics_setup.json"
    };
    for (const auto& p : candidates) {
        if (std::filesystem::exists(p)) {
            return p;
        }
    }
    return user_hint;
}

/**
 * @brief 从 docs/physics_setup.json 动态加载 CHANNEL_CONFIGS 字典
 * 彻底告别在 C++ 中硬编码信道映射，保证 docs/ 为单一可信源！
 */
[[nodiscard]] inline std::map<std::string, std::vector<MesonChannelMapping>>
load_channel_configs_from_docs(const std::filesystem::path& json_hint = "docs/physics_setup.json") {
    const auto resolved = resolve_physics_setup_path(json_hint);
    const std::string content = iodata::read_file_to_string(resolved);
    if (content.empty()) {
        throw std::runtime_error("无法读取物理配置文件: " + resolved.string());
    }

    std::map<std::string, std::vector<MesonChannelMapping>> channel_defs;

    // 1. 定位 "CHANNEL_CONFIGS"
    const size_t cfg_pos = content.find("\"CHANNEL_CONFIGS\"");
    if (cfg_pos == std::string::npos) {
        throw std::runtime_error("配置文件中未找到 'CHANNEL_CONFIGS': " + resolved.string());
    }

    const size_t obj_start = content.find('{', cfg_pos);
    if (obj_start == std::string::npos) {
        return channel_defs;
    }

    // 2. 查找匹配的闭合大括号，限定作用域
    size_t pos = obj_start + 1;
    int brace_count = 1;
    size_t obj_end = std::string::npos;
    while (pos < content.size()) {
        if (content[pos] == '{') brace_count++;
        else if (content[pos] == '}') {
            brace_count--;
            if (brace_count == 0) {
                obj_end = pos;
                break;
            }
        }
        pos++;
    }

    if (obj_end == std::string::npos) {
        throw std::runtime_error("CHANNEL_CONFIGS JSON 结构不完整");
    }

    const std::string_view cfg_scope(content.data() + obj_start, obj_end - obj_start + 1);

    // 3. 在 cfg_scope 内解析形如 "AV": [ ... ]
    size_t cursor = 0;
    while (cursor < cfg_scope.size()) {
        // 查找下一个键："channel_name": [
        size_t key_quote_start = cfg_scope.find('"', cursor);
        if (key_quote_start == std::string_view::npos) break;
        size_t key_quote_end = cfg_scope.find('"', key_quote_start + 1);
        if (key_quote_end == std::string_view::npos) break;

        std::string channel_name(cfg_scope.substr(key_quote_start + 1, key_quote_end - key_quote_start - 1));

        // 查找紧随其后的 '[' 数组开始
        size_t arr_start = cfg_scope.find('[', key_quote_end + 1);
        if (arr_start == std::string_view::npos) break;

        // 查找对应的 ']' 数组结束
        size_t arr_end = cfg_scope.find(']', arr_start + 1);
        if (arr_end == std::string_view::npos) break;

        std::string_view arr_scope = cfg_scope.substr(arr_start + 1, arr_end - arr_start - 1);

        // 在 arr_scope 内解析所有的 { "type": "...", "dir": "..." }
        std::vector<MesonChannelMapping> mappings;
        size_t item_pos = 0;
        while (item_pos < arr_scope.size()) {
            size_t item_start = arr_scope.find('{', item_pos);
            if (item_start == std::string_view::npos) break;
            size_t item_end = arr_scope.find('}', item_start);
            if (item_end == std::string_view::npos) break;

            std::string_view item_str = arr_scope.substr(item_start + 1, item_end - item_start - 1);

            // 提取 "type"
            std::string type_val;
            size_t type_tag = item_str.find("\"type\"");
            if (type_tag != std::string_view::npos) {
                size_t q1 = item_str.find('"', type_tag + 6);
                if (q1 != std::string_view::npos) {
                    size_t q2 = item_str.find('"', q1 + 1);
                    if (q2 != std::string_view::npos) {
                        type_val = std::string(item_str.substr(q1 + 1, q2 - q1 - 1));
                    }
                }
            }

            // 提取 "dir"
            std::string dir_val;
            size_t dir_tag = item_str.find("\"dir\"");
            if (dir_tag != std::string_view::npos) {
                size_t q1 = item_str.find('"', dir_tag + 5);
                if (q1 != std::string_view::npos) {
                    size_t q2 = item_str.find('"', q1 + 1);
                    if (q2 != std::string_view::npos) {
                        dir_val = std::string(item_str.substr(q1 + 1, q2 - q1 - 1));
                    }
                }
            }

            if (!type_val.empty() && !dir_val.empty()) {
                mappings.push_back(MesonChannelMapping{type_val, dir_val});
            }

            item_pos = item_end + 1;
        }

        channel_defs[channel_name] = std::move(mappings);
        cursor = arr_end + 1;
    }

    return channel_defs;
}

/**
 * @brief 从 docs/physics_setup.json 动态加载 CONDENSATE_CONFIGS 列表
 */
[[nodiscard]] inline std::vector<CondensateConfig>
load_condensate_configs_from_docs(const std::filesystem::path& json_hint = "docs/physics_setup.json") {
    const auto resolved = resolve_physics_setup_path(json_hint);
    const std::string content = iodata::read_file_to_string(resolved);
    if (content.empty()) {
        throw std::runtime_error("无法读取物理配置文件: " + resolved.string());
    }

    std::vector<CondensateConfig> configs;

    const size_t cfg_pos = content.find("\"CONDENSATE_CONFIGS\"");
    if (cfg_pos == std::string::npos) {
        return configs;
    }

    const size_t arr_start = content.find('[', cfg_pos);
    const size_t arr_end = content.find(']', arr_start);
    if (arr_start == std::string::npos || arr_end == std::string::npos) {
        return configs;
    }

    const std::string_view arr_str(content.data() + arr_start, arr_end - arr_start + 1);
    size_t item_pos = 0;
    while (item_pos < arr_str.size()) {
        size_t obj_s = arr_str.find('{', item_pos);
        if (obj_s == std::string_view::npos) break;
        size_t obj_e = arr_str.find('}', obj_s);
        if (obj_e == std::string_view::npos) break;

        std::string_view obj_str = arr_str.substr(obj_s, obj_e - obj_s + 1);

        auto extract_num = [&](std::string_view key, double& out_val) {
            size_t kpos = obj_str.find(key);
            if (kpos == std::string_view::npos) return;
            size_t colon = obj_str.find(':', kpos);
            if (colon == std::string_view::npos) return;
            size_t val_start = colon + 1;
            while (val_start < obj_str.size() && (obj_str[val_start] == ' ' || obj_str[val_start] == '\t')) val_start++;
            size_t val_end = obj_str.find_first_of(",}\r\n ", val_start);
            if (val_end == std::string_view::npos) val_end = obj_str.size();
            iodata::parse_double(obj_str.substr(val_start, val_end - val_start), out_val);
        };

        CondensateConfig c;
        extract_num("\"beta\"", c.beta);
        extract_num("\"ml\"", c.ml);
        extract_num("\"ms\"", c.ms);
        extract_num("\"mres\"", c.mres);
        extract_num("\"pbp_l\"", c.pbp_l);
        extract_num("\"pbp_l_err\"", c.pbp_l_err);
        extract_num("\"pbp_s\"", c.pbp_s);
        extract_num("\"pbp_s_err\"", c.pbp_s_err);
        extract_num("\"zm\"", c.zm);
        extract_num("\"pbp_rm\"", c.pbp_rm);
        extract_num("\"pbp_rm_err\"", c.pbp_rm_err);

        configs.push_back(c);
        item_pos = obj_e + 1;
    }

    return configs;
}

} // namespace lqcd
