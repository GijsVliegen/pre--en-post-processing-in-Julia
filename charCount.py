import csv

# Script om gemiddelde aantal lijnen te tellen in een blok Python of Julia code.
result = open("results.csv", 'w', newline='')
writer = csv.writer(result)
code = open("input.txt", 'r',encoding='utf-8')

writer.writerow(["PYTHON"])
writer.writerow([])
writer.writerow(["number of lines", "number of characters", "avg chars per line", "number of brackets", "average brackets per line"])



assert code.readline() == "PYTHON\n"

while True:
    line = code.readline()
    if (line.startswith("FUNCTION")):
        lengths = []
        brackets = 0
        nb_brackets = 0
        while True:

            line = code.readline().replace("\n", '')
            if line.startswith("END"):
                writer.writerow([
                    len(lengths), # number of lines
                    sum(lengths), # number of characters
                    sum(lengths)/len(lengths), # avg chars per line
                    nb_brackets, # number of brackets
                    nb_brackets/len(lengths) # avg pairs of brackets per line
                ])
                break

            # remove whitespace
            line = line.replace(' ','')

            # handle line breaks and brackets
            for i in range(0, len(line)):
                if line[i] == "(" or line[i] == "{" or line[i] == "[":
                    nb_brackets += 1
                if line[i] == "(":
                    brackets += 1
                elif line[i] == ")":
                    brackets -= 1
            if brackets == 0:
                lengths.append(len(line))
            elif brackets > 0:
                if len(lengths) == 0:
                    lengths.append(0)
                lengths[-1] += len(line)

    if line.startswith("JULIA"):
        break

writer.writerow(["JULIA"])
writer.writerow([])
writer.writerow(["number of lines", "number of characters", "avg chars per line", "number of brackets", "average brackets per line"])

# Process Julia functions
while True:
    line = code.readline()

    if line.startswith("FUNCTION"):
        lengths = []
        brackets = 0
        nb_brackets = 0
        while True:
            line = code.readline().replace("\n", '')
            if line.startswith("END"):
                writer.writerow([
                    len(lengths), # number of lines
                    sum(lengths), # number of characters
                    sum(lengths)/len(lengths), # avg chars per line
                    nb_brackets, # number of brackets
                    nb_brackets/len(lengths) # avg pairs of brackets per line
                ])
                break

            # remove whitespace
            line = line.replace(' ', '')

            lengths.append(len(line))
            if "{" in line or "(" in line or "[" in line:
                nb_brackets += 1

    if line.startswith("FILE END"):
        break


result.close()