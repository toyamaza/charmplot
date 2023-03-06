#!/usr/bin/env python
import yaml


def open_yaml(file):
    # Open the YAML file
    with open(file, "r") as file:
        # Load the YAML document into a Python object
        data = yaml.safe_load(file)
        # Use the data as needed
        return data


def add_zero(string):
    if abs(float(string)) < 1e-5:
        return "0.0"
    n = 0
    for c in string:
        if c != "0" and c != "." and c != "-":
            n += 1
    if n < 2:
        string += "0"

    n_zeros = 0
    for c in string[::-1]:
        if c == "0":
            n_zeros += 1
        else:
            break
    if n == 1:
        n_zeros -= 1
    for _ in range(n_zeros):
        string = string[:-1]

    return string


def main():

    for meson in ["dplus", "dstar"]:
        dplus_minus = open_yaml(f"hepdata/{meson}_minus.yaml")
        dplus_plus = open_yaml(f"hepdata/{meson}_plus.yaml")

        fout = open(f"{meson}.txt", "w")
        for index, x in enumerate(dplus_minus["dependent_variables"][1]["values"][0]["errors"]):
            # NP name
            label = x["label"]
            if label == "stat" or "Uncor" in label:
                continue

            string = x["label"]
            if string.startswith("gamma_stat_"):
                string = string.replace("gamma_stat_", "$\\gamma$(")
                string += ")"
            string = string.replace("_", "\\_")

            # minus
            for i in range(5):
                y = dplus_minus["dependent_variables"][1]["values"][i]["errors"][index]
                minus = f"{y['asymerror']['minus']:f}"
                plus = f"{y['asymerror']['plus']:f}"

                # add trailing zero
                minus = add_zero(minus)
                plus = add_zero(plus)

                # add sign
                if minus[0] != "-":
                    minus = "+" + minus
                if plus[0] != "-":
                    plus = "+" + plus
                string += f" & $^{{{plus}}}_{{{minus}}}$"

            # plus
            for i in range(5):
                found = False
                for y in dplus_plus["dependent_variables"][1]["values"][i]["errors"]:
                    if y["label"] == label:
                        found = True
                        break
                assert found, f"missing NP {label}"
                minus = f"{y['asymerror']['minus']:f}"
                plus = f"{y['asymerror']['plus']:f}"

                # add trailing zero
                minus = add_zero(minus)
                plus = add_zero(plus)

                # add sign
                if minus[0] != "-":
                    minus = "+" + minus
                if plus[0] != "-":
                    plus = "+" + plus
                string += f" & $^{{{plus}}}_{{{minus}}}$"

            string += " \\\\ \\hline"
            print(string)
            fout.write(string + "\n")

        for label in [f"Uncorr_{meson}_mu_W{charge}_{b}" for charge in ["minus", "plus"] for b in range(1, 6)]:
            string = label.replace("_", "\\_")

            for i in range(5):
                found = False
                for y in dplus_minus["dependent_variables"][1]["values"][i]["errors"]:
                    if y["label"] == label:
                        found = True
                        break
                if found:
                    minus = f"{y['asymerror']['minus']:f}"
                    plus = f"{y['asymerror']['plus']:f}"
                else:
                    minus = f"{0:f}"
                    plus = f"{0:f}"

                # add trailing zero
                minus = add_zero(minus)
                plus = add_zero(plus)

                # add sign
                if minus[0] != "-":
                    minus = "+" + minus
                if plus[0] != "-":
                    plus = "+" + plus
                string += f" & $^{{{plus}}}_{{{minus}}}$"

            for i in range(5):
                found = False
                for y in dplus_plus["dependent_variables"][1]["values"][i]["errors"]:
                    if y["label"] == label:
                        found = True
                        break
                if found:
                    minus = f"{y['asymerror']['minus']:f}"
                    plus = f"{y['asymerror']['plus']:f}"
                else:
                    minus = f"{0:f}"
                    plus = f"{0:f}"

                # add trailing zero
                minus = add_zero(minus)
                plus = add_zero(plus)

                # add sign
                if minus[0] != "-":
                    minus = "+" + minus
                if plus[0] != "-":
                    plus = "+" + plus
                string += f" & $^{{{plus}}}_{{{minus}}}$"

            string += " \\\\ \\hline"
            print(string)
            fout.write(string + "\n")

        # close file
        fout.close


if __name__ == "__main__":

    # run
    main()
