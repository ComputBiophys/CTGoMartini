import os
import subprocess

def run_tests():
    data_path = os.path.split(__file__)[0]
    try:
        subprocess.run(["pytest", data_path], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Test failed with return code {e.returncode}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    run_tests()
