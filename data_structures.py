#Create a list containing the of the names of your fellow trainees
participant_list = ['Andrés', 'Hongyang',
                    'Phil', 'Jen', 'Harry',
                    'Maisha', 'Sophie', 'Hebe']

#• Add the names of today's trainers to the list using append
participant_list.append("David")
participant_list.append("Charlie")
participant_list.append("Kevin")
print(participant_list)

#• Select the 3rd and 5th names from your list
print(participant_list[2])
print(participant_list[4])

#• Sort your list & select the 3rd to the 5th names from your list
participant_list.sort()
print(participant_list)
print(participant_list[2:5])

#Iterate over the list and set the course participants names to keys in a dictionary with the value as 'participant' for each of them
participant_dictionary = {}
for participant in participant_list:
    participant_dictionary[participant] = "participant"
    


#Use a for loop to iterate over your dictionary and print the values of the keys only if they are trainers
trainer_list = ['Kevin', 'Charlie', 'David']



print(participant_dictionary)

for participant in participant_dictionary.values():
    print(participant)
    
for participant in participant_dictionary.items():
    print(participant)
    

for participant in participant_dictionary.keys():
    if participant in trainer_list:
        print(participant)
    else:
        continue


#For select the first 2 letters of the string in the third value of your list
print(participant_list[2][0:2])