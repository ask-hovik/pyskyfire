import os

def count_nonempty_lines_in_src(src_path="src"):
    total_lines = 0
    for root, _, files in os.walk(src_path):
        for file in files:
            if file.endswith(".py"):  # Only count Python files
                file_path = os.path.join(root, file)
                with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
                    for line in f:
                        if line.strip():  # Non-empty line
                            total_lines += 1
    return total_lines

if __name__ == "__main__":
    count = count_nonempty_lines_in_src("src")
    print(f"Total non-empty lines of code in src/: {count}")
