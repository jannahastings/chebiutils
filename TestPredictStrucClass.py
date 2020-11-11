import os
os.chdir('/Users/hastingj/Work/Python/chemont/chemont-struc/experiments')

import pandas
import matplotlib.pyplot as plt
from sklearn import model_selection
from sklearn.metrics import classification_report
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
import seaborn as sns
from matplotlib_venn import venn3, venn3_circles
import pickle
from collections import Counter
import numpy as np


from chebidblite.learnhelper import ChebiDataPreparer
from chebidblite.learnhelper import ChebiOntologySubsetter

import csv


RANDOM_STATE = 42

# ChebiDataPreparer from dblite/learnhelper
dprep = ChebiDataPreparer()

# Load the dataset
problem_sizes = ['25x25',
                 '25x100',
                 '25x200',
                 '25x500',
                 '100x25',
                 '100x100',
                 '100x200',
                 '100x500',
                 '150x25',
                 '150x100',
                 '150x500',
                 '200x25',
                 '200x100',
                 '200x200',
                 '200x500',
                 '500x25',
                 '500x100',
                 '500x200']

for problem_size in problem_sizes:
    print("#### PROCESSING problem size",problem_size,"####")
    
    [n_classes,n_members] = problem_size.split('x')
    
    datafilename = 'data/chemdata'+problem_size+'x1024.pkl'
    if not os.path.exists(datafilename):
        # Uncomment this to generate a new dataset (slow)
        #chemdata = dprep.getDataForLearning(class_size=int(n_members),
        #                                    number_to_select=int(n_classes))
        #with open(datafilename,'wb') as output:
        #    pickle.dump(chemdata,output)
        print("File ",datafilename," not found. Skipping...")
        continue
    else:
        chemdata = pickle.load(open(datafilename,'rb'))
    
    # Verify the number of classes is as you would expect
    if not len(set(chemdata.iloc[:,[0]][0].values)) == int(n_classes):
        print("NOTE: MISMATCH IN NUMBER OF CLASSES: ",len(set(chemdata.iloc[:,[0]][0].values)),"DOES NOT MATCH EXPECTED",n_classes)
    
    X = chemdata.iloc[:,3:1027] # Get just the fingerprint columns
    y = chemdata.iloc[:,[0]][0].values # Get just the class ID column
    
    # Get a training and testing split of the data (20% for testing)
    X_train, X_test, y_train, y_test = model_selection.train_test_split(X, y, test_size=0.20, random_state=RANDOM_STATE)
    
    # COMPARE FLAT, SHALLOW ALGORITHMS
    models = []
    models.append(('LR', LogisticRegression(solver='liblinear', multi_class='ovr')))
    models.append(('KNN', KNeighborsClassifier()))
    models.append(('CART', DecisionTreeClassifier()))
    models.append(('RF', RandomForestClassifier()))
    models.append(('NB', GaussianNB()))
    models.append(('SVM-linear', SVC(kernel='linear', gamma='auto')))
    models.append(('SVM-rbf', SVC(kernel='rbf', gamma='auto')))
    models.append(('SVM-sigmoid', SVC(kernel='sigmoid', gamma='auto')))
    models.append(('LDA', LinearDiscriminantAnalysis()))
    
    classif_reports = {}
    # Make predictions on validation datasets
    for name, model in models:
        model.fit(X_train, y_train)
        
        predictions = model.predict(X_test)
        
        print("Accuracy of",name, ":", accuracy_score(y_test, predictions))
        print("Weighted average precision of",name, ":", precision_score(y_test, predictions, average='weighted'))
        print("Weighted average recall of",name, ":", recall_score(y_test, predictions, average='weighted'))
        
        classifres = classification_report(y_test, predictions, output_dict=True)
        classif_reports[name] = classifres
    
    # Save the results
    with open('data/classif_reports_'+problem_size+'x1024.pkl','wb') as output:
        pickle.dump(classif_reports,output)
       


# Overall analyses of the results, depends on parsing the result files

# Load the average scores per classifier type and per problem size

overall_results = {}
overall_results['problem_descr'] = []
overall_results['classes'] = []
overall_results['members'] = []
overall_results['algorithm'] = []
overall_results['av_accuracy'] = []
overall_results['av_precision'] = []
overall_results['av_recall'] = []
overall_results['av_f1'] = []

for problem_size in problem_sizes:
    [n_classes,n_members] = problem_size.split('x')
    
    resultsfilename = 'data/classif_reports_'+problem_size+'x1024.pkl'
    if not os.path.exists(resultsfilename):
        print("### CLASSIF REPORT NOT FOUND FOR PROBLEM SIZE ",problem_size)
        continue

    with open(resultsfilename,'rb') as input:
        classif_reports = pickle.load(input)
    
    groups = [k for k in classif_reports.keys()]
    # Sanity check size of the report
    chebi_ids = [k for k in classif_reports[groups[1]].keys() if "CHEBI" in k]
    
    if not len(set(chebi_ids)) == int(n_classes):
        print("### CLASSIF SIZE MISMATCH: ",len(set(chebi_ids))," vs expected ",n_classes)
    
    # Compare average performance metrics 
    for model in classif_reports:
        accuracy = classif_reports[model]['accuracy']
        av_prec = classif_reports[model]['weighted avg']['precision']
        av_recall = classif_reports[model]['weighted avg']['recall']
        av_f1 = classif_reports[model]['weighted avg']['f1-score']
        overall_results['problem_descr'].append(problem_size)
        overall_results['classes'].append(int(n_classes))
        overall_results['members'].append(int(n_members))
        overall_results['algorithm'].append(model)
        overall_results['av_accuracy'].append(accuracy)
        overall_results['av_precision'].append(av_prec)
        overall_results['av_recall'].append(av_recall)
        overall_results['av_f1'].append(av_f1)
            
        
# Add the LSTM problem size scores 
for problem_size in ["100x100","100x500","500x100"]:
    [n_classes,n_members] = problem_size.split('x')
    with (open("results/lstm-results-classes-"+problem_size+".csv", 'r')) as csvfile:
        csvreader = csv.DictReader(csvfile)
        f1_scores = []
        for row in csvreader:
            f1_score = float(row['f1'])
            f1_scores.append(f1_score)
        av_f1 = np.mean(f1_scores)
        overall_results['problem_descr'].append(problem_size)
        overall_results['classes'].append(int(n_classes))
        overall_results['members'].append(int(n_members))
        overall_results['algorithm'].append('LSTM')
        overall_results['av_accuracy'].append(av_f1)
        overall_results['av_precision'].append(av_f1)
        overall_results['av_recall'].append(av_f1)
        overall_results['av_f1'].append(av_f1)
        
            
df = pandas.DataFrame.from_dict(overall_results)    
df = df.query("algorithm != 'SVM'")
    

df_lineplot = df.query("members == 25 | members == 100 | members == 200")
# Now plot the averages of all the F1 scores by size parameters
ax = sns.lineplot(x="classes", y="av_f1", 
                      hue="algorithm",
                      style="members",
                      data=df )
plt.ylim(0.1, 1.1)
plt.xlim(0,700)
plt.legend(loc='upper right')
plt.setp(ax.get_legend().get_texts(), fontsize='7') # for legend text




# Now look at the different sub-classes that are selected and predict for them
    
f1_results = {}
f1_results['problem_descr'] = []
f1_results['classes'] = []
f1_results['members'] = []
f1_results['algorithm'] = []
f1_results['class_id'] = []
f1_results['f1_score'] = []

for problem_size in problem_sizes:
    [n_classes,n_members] = problem_size.split('x')
    
    resultsfilename = 'data/classif_reports_'+problem_size+'x1024.pkl'
    if not os.path.exists(resultsfilename):
        print("### CLASSIF REPORT NOT FOUND FOR PROBLEM SIZE ",problem_size)
        continue

    with open(resultsfilename,'rb') as input:
        classif_reports = pickle.load(input)
    
    for model in classif_reports: 
        classifres = classif_reports[model]
        chebi_ids = [k for k in classifres.keys() if "CHEBI" in k]
        f1s = [classifres[k]['f1-score'] for k in classifres.keys() if "CHEBI" in k]
        
        for (chebi_id,f1) in zip(chebi_ids,f1s):
            f1_results['problem_descr'].append(problem_size)
            f1_results['classes'].append(int(n_classes))
            f1_results['members'].append(int(n_members))
            f1_results['algorithm'].append(model)
            f1_results['class_id'].append(chebi_id)
            f1_results['f1_score'].append(f1)
            
    

# Add the LSTM problem size scores 

for problem_size in ["100x100","100x500","500x100"]:
    [n_classes,n_members] = problem_size.split('x')
    with (open("results/lstm-results-classes-"+problem_size+".csv", 'r')) as csvfile:
        csvreader = csv.DictReader(csvfile)
        for row in csvreader:
            class_id = row['id']
            f1_score = float(row['f1'])
            
            f1_results['problem_descr'].append('100x500')
            f1_results['classes'].append(int(n_classes))
            f1_results['members'].append(int(n_members))
            f1_results['algorithm'].append('LSTM')
            f1_results['class_id'].append(class_id)
            f1_results['f1_score'].append(f1_score)
        
dff1 = pandas.DataFrame.from_dict(f1_results)    


# Draw a box plot of the F1 scores for the different algorithms
dff1_subs = dff1.query('classes==25 & members==100')
b = sns.boxplot(x='algorithm', y='f1_score', data=dff1_subs)
b.set_title("F1 scores 25x100")
b.tick_params(labelsize=7)

dff1_subs = dff1.query('classes==100 & members==100')
b = sns.boxplot(x='algorithm', y='f1_score', data=dff1_subs)
b.set_title("F1 scores 100x100")
b.tick_params(labelsize=7)

dff1_subs = dff1.query('classes==100 & members==500')
b = sns.violinplot(x='algorithm', y='f1_score', data=dff1_subs)
b.set_title("F1 scores 100x500")
b.tick_params(labelsize=7)


dff1_subs = dff1.query('classes==500 & members==25')
b = sns.boxplot(x='algorithm', y='f1_score', data=dff1_subs)
b.set_title("F1 scores 500x25")
b.tick_params(labelsize=7)

dff1_subs = dff1.query('classes==500 & members==100')
b = sns.boxplot(x='algorithm', y='f1_score', data=dff1_subs)
b.set_title("F1 scores 500x100")
b.tick_params(labelsize=7)


dff1_subs = dff1.query('classes==500 & members==100')
b = sns.violinplot(x='algorithm', y='f1_score', data=dff1_subs)
b.set_title("F1 scores 500x100")
b.tick_params(labelsize=7)

dff1_subs = dff1.query('classes==25 & ( members==100 | members == 500 )')
ax = sns.violinplot(x='algorithm', y="f1_score", hue="members",
                    data=dff1_subs, palette="muted", split=True)
ax.set_title("F1 scores 25x100 and 25x500")
ax.tick_params(labelsize=7)

dff1_subs = dff1.query('classes==100 & ( members==100 | members ==500 )')
ax = sns.violinplot(x='algorithm', y="f1_score", hue="members",
                    data=dff1_subs, palette="muted", split=True)
ax.set_title("F1 scores 100x100 and 100x500")
ax.tick_params(labelsize=7)

dff1_subs = dff1.query('classes==200 & ( members==100 | members ==500 )')
ax = sns.violinplot(x='algorithm', y="f1_score", hue="members",
                    data=dff1_subs, palette="muted", split=True)
ax.set_title("F1 scores 200x100 and 200x500")
ax.tick_params(labelsize=7)


# Taking the fraction with F1 > 0.8 in each classifier. 
# Do they contain the same ChEBI classes, or different? 
# (I.e. can we achieve better overall performance by combining)
# (classifiers for different classes? )
for n_classes in [100,200,500]:
    for n_members in [25,100,200]:
        dff1_subs = dff1.query('classes=='+str(n_classes)+' &  members=='+str(n_members)+' & f1_score > 0.8')
    
        class_lists = dff1_subs.groupby('algorithm')['class_id'].apply(list)
    
        # Take the best 3 and draw a VENN: LR, LDA, RF (SVM)
        # Pick from the others for another VENN: NB, CART, KNN
        
        classifiers_for_venn = ['LR','LDA','RF'] # choose any 3 (LSTM has restricted problem sizes available)
        
        data_for_venn = [ set(class_lists[c]) for c in classifiers_for_venn ]
        
        vd3=venn3(data_for_venn,
                  set_labels=(classifiers_for_venn),
                  set_colors=('#c4e6ff', '#F4ACB7','#9D8189'), 
                  alpha = 0.8)
        venn3_circles(data_for_venn, linestyle='-.', linewidth=2, color='grey')
        for text in vd3.set_labels:
            text.set_fontsize(16);
        plt.title('Top-performing classes'+str(n_classes)+'x'+str(n_members),fontname='Times New Roman',fontweight='bold',fontsize=20,
                  pad=30,backgroundcolor='#cbe7e3',color='black',style='italic');
        plt.show()
        

### Now get the top-ranked classes (F1>0.8) for all approaches

best_classes_counts = {}
best_classes_counts["problem_size"] = []
best_classes_counts["classes"] = []
best_classes_counts["members"] = []
best_classes_counts["chebi_id"] = []
best_classes_counts["count"] = []
for problem_size in problem_sizes:
    [n_classes,n_members] = problem_size.split('x')
    dff1_subs = dff1.query('classes=='+str(n_classes)+' &  members=='+str(n_members)+' & f1_score > 0.8')

    class_lists = dff1_subs.groupby('algorithm')['class_id'].apply(list)

    print(problem_size, class_lists)
    
    class_lists_in = [ i for c in class_lists.keys() for i in class_lists[c] ]
    
    ccc = Counter(class_lists_in).keys() 
    ccc_values = Counter(class_lists_in).values() # counts the elements' frequency
    print(ccc, ccc_values)
    
    for cc,vv in zip(ccc,ccc_values): 
        best_classes_counts["problem_size"].append(problem_size)
        best_classes_counts["classes"].append(int(n_classes))
        best_classes_counts["members"].append(int(n_members))
        best_classes_counts["chebi_id"].append(cc)
        best_classes_counts["count"].append(vv)

dff2 = pandas.DataFrame.from_dict(best_classes_counts)    

for n_classes in [100,200,500]:
    for n_members in [25,100,200]:
        dff2_subs = dff2.query('classes=='+str(n_classes)+' &  members=='+str(n_members))

        data_ag = dff2_subs.groupby('chebi_id')['count'].sum()
        data_ag = data_ag.sort_values(ascending=False)

        subsetter = ChebiOntologySubsetter()
        subsetter.createSubsetFor(data_ag.keys())
        colour_numbers = data_ag/max(data_ag)
        subsetter.printImageOfSubset(image_name="results/network_best_classes"+str(n_classes)+"x"+str(n_members)+".png", colour_nums=colour_numbers)




### Now get the mean F1 score for all classes across all approaches
dff3 = dff1.groupby('class_id')['f1_score'].mean()
dff3 = dff3.sort_values(ascending=False)
class_ids = dff3.index.values
dff3 = pandas.DataFrame(dff3)
dff3['class_name'] = [e.chebi_name for e in [dprep.db.getEntity(c) for c in class_ids] ]
dff3['class_members'] = [len(dprep.classes_to_members[c]) for c in class_ids]

# Include tables in the 
print(dff3.head(20).to_latex())
print(dff3.tail(20).to_latex())


# Just the LSTMs?
dff4 = dff1.query("algorithm == 'LSTM' ").groupby('class_id')['f1_score'].mean()
dff4 = dff4.sort_values(ascending=False)
class_ids = dff4.index.values
dff4 = pandas.DataFrame(dff4)
dff4['class_name'] = [e.chebi_name for e in [dprep.db.getEntity(c) for c in class_ids] ]
dff4['class_members'] = [len(dprep.classes_to_members[c]) for c in class_ids]

# Include tables in the paper
print(dff4.head(20).to_latex())
print(dff4.tail(20).to_latex())



# Visualisation across the full ontology?


subsetter = ChebiOntologySubsetter()
subsetter.createSubsetFor(dff4.index.values)
colour_numbers = dff4['f1_score']
subsetter.printImageOfSubset(image_name="results/lstm_best_classes.png", colour_nums=colour_numbers)


# Display the F1 scores as colours on a network graph representation of the ChEBI ontology.

#subsetter = ChebiOntologySubsetter()
#subsetter.createSubsetFor(lr_200x100x1024_chebis)
#colour_numbers = {k:v for (k,v) in zip(lr_200x100x1024_chebis,lr_200x100x1024_f1s)}
#subsetter.printImageOfSubset(image_name="results/network_lr_200x100x1024.png", colour_nums=colour_numbers)








#### Number of classes for N numbers of members: 
    

sizes = [10,25,50,100,200,500]
counts = []
for size in sizes:
    count =  len(  [c for c in dprep.classes_to_members.keys() if len(dprep.classes_to_members[c])>size] )
    counts.append(count)
y_pos = np.arange(len(sizes)) 
# Create bars
plt.bar(y_pos, counts)
# Create names on the x-axis
plt.xticks(y_pos, sizes)
plt.xlabel("Number of members (M)")
plt.ylabel("Number of classes")
# Show graphic
plt.show()

    
### Make a plot of number of members (actually in ChEBI) vs. F1 score
### Do smaller (more specific) classes perform better on average?



# Density Plot and Histogram of F1 scores
        #sns.distplot(f1s, hist=True, kde=True,
        #            bins=int(180/5), color = 'darkblue',
        #            hist_kws={'edgecolor':'black'},
        #            kde_kws={'linewidth': 4})
        #plt.title('Density plot of F1 scores for '+model+","+problem_size)
        #plt.ylabel('Scores')
        #plt.show()
