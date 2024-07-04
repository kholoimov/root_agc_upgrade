import json

data = json.load(open('my_json.json'))

file = open("myjson.json", "w")
file.write(json.dumps(data, indent=4))
file.close()