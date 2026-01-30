def args_to_str(args_dict, pair_seperator="_"):
    parts = []
    for k, v in args_dict.items():
        parts.append(f"{k}={v}")
    return pair_seperator.join(parts)
