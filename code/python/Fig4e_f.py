#!/usr/bin/env python
# coding: utf-8

# In[32]:


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import random
plt.rcParams.update({
    "font.size": 10,              # base font
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 12
})
# Read the data
Topo_Network_5 = pd.read_csv('./data/RACIPE_practice/TS_2.topo', sep=r'\s+')

def Topo_adj(dp):
    n = len(set(dp['Source']))
    adj = pd.DataFrame([[0]*n]*n)
    adj.columns = set(dp.loc[:, 'Source'])
    adj.index = set(dp.loc[:, 'Source'])
    for i in range(len(dp)):
        if dp.iloc[i, 2] == 1:
            adj.loc[dp.iloc[i, 0], dp.iloc[i, 1]] = 1
        elif dp.iloc[i, 2] == 2:
            adj.loc[dp.iloc[i, 0], dp.iloc[i, 1]] = -1
    return adj


g1 = sns.clustermap(
    Topo_adj(Topo_Network_5),
    annot=False,
    cmap="coolwarm",
    linecolor='black',
    linewidths=1,
    figsize=(4, 3)
    
)

g1.fig.suptitle(
    "Adjacency Matrix",
    fontsize=12,
    y=1.02
)
plt.setp(g1.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
plt.setp(g1.ax_heatmap.get_yticklabels())
g1.fig.savefig(
    "./figures/Fig4/Fig4Ei_adjacency.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close(g1.fig)


# Influence Calculation
def influence(dp, j):
    n = len(set(dp['Source']))
    infl = pd.DataFrame([[0]*n]*n)
    infl.columns = set(dp.loc[:, 'Source'])
    infl.index = set(dp.loc[:, 'Source'])
    adj = Topo_adj(dp)

    maxi = abs(adj)
    pr = adj
    mr = abs(adj)
    infl = adj
    for i in range(1, j):
        pr = adj @ pr
        mr = maxi @ mr
        qwe = pr / mr

        infl = infl + qwe.replace(np.nan, 0)

    infl = infl / j
    infl = infl.round(2)
    return infl





tup = influence(Topo_Network_5, 8).sort_values(by='mir200', ascending=True)
tupi = tup.T.sort_values(by='mir200', ascending=True)
print(tupi)

g2 = sns.clustermap(
    influence(Topo_Network_5, 8),
    annot=False,
    cmap='seismic',
    vmax=1,
    vmin=-1,
    linecolor='black',
    linewidths=1,
    figsize=(4, 3)
)

g2.fig.suptitle(
    "Influence Matrix",
    fontsize=12,
    y=1.02
)
plt.setp(g2.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
plt.setp(g2.ax_heatmap.get_yticklabels())
g2.fig.savefig(
    "./figures/Fig4/Fig4Eii_influence.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close(g2.fig)
# Identifying Teams
Team1 = []
Team2 = []
infl = influence(Topo_Network_5, 8)
for i in infl.index:
    if infl.loc[i, 'mir200'] >= 0:
        Team1.append(i)
    elif infl.loc[i, 'mir200'] <= 0:
        Team2.append(i)

print(Team1)
print(Team2)
# Calculating Team Scores
TeamScores = pd.DataFrame(0.0, columns=[0, 1], index=[0, 1])
teams = [Team1, Team2]
for t1 in range(len(teams)):
    for t2 in range(len(teams)):
        Tsc = 0
        count = 0
        for i in teams[t1]:
            for j in teams[t2]:
                Tsc += tupi.loc[i, j]
                count += 1
        TeamScores.loc[t1, t2] = Tsc / count

TeamScores.columns = ['TeamB', 'TeamA']
TeamScores.index = ['TeamB', 'TeamA']
print(TeamScores)
teamscore_wt = abs(TeamScores).sum(axis=1).sum(axis=0) / 4
print(teamscore_wt)
# More Randomization (1000 times)
log = []
dp = Topo_Network_5

for oi in range(1000):
    edges = list(dp['Type'])
    random.shuffle(edges)
    dp['Type'] = edges

    n = len(set(dp['Source']))
    adj = pd.DataFrame([[0.0]*n]*n)
    adj.columns = set(dp.loc[:, 'Source'])
    adj.index = set(dp.loc[:, 'Source'])

    for i in range(len(dp)):
        if dp.iloc[i, 2] == 1:
            adj.loc[dp.iloc[i, 0], dp.iloc[i, 1]] = 1
        elif dp.iloc[i, 2] == 2:
            adj.loc[dp.iloc[i, 0], dp.iloc[i, 1]] = -1

    n = len(set(dp['Source']))
    infl = pd.DataFrame([[0.0]*n]*n)
    infl.columns = set(dp.loc[:, 'Source'])
    infl.index = set(dp.loc[:, 'Source'])

    maxi = abs(adj)
    pr = adj
    mr = abs(adj)
    infl = adj
    for i in range(1, 8):
        pr = adj @ pr
        mr = maxi @ mr
        qwe = pr / mr

        infl = infl + qwe.replace(np.nan, 0)

    infl = infl / 8
    infl

    Team1 = []
    Team2 = []
    for i in infl.index:
        if infl.loc[i, 'mir200'] >= 0:
            Team1.append(i)
        elif infl.loc[i, 'mir200'] <= 0:
            Team2.append(i)

    TeamScores = pd.DataFrame(0.0, columns=[0, 1], index=[0, 1])
    
    teams = [Team1, Team2]
    for t1 in range(len(teams)):
        for t2 in range(len(teams)):
            Tsc = 0
            count = 1
            for i in teams[t1]:
                for j in teams[t2]:
                    Tsc += infl.loc[i, j]
                    count += 1
            TeamScores.loc[t1, t2] = Tsc / count

    TeamScores.columns = ['TeamB', 'TeamA']
    TeamScores.index = ['TeamB', 'TeamA']
    TeamScores

    log.append(abs(TeamScores).sum(axis=1).sum(axis=0) / 4)

qwer = pd.DataFrame(log)
qwer[1] = qwer.index
print(qwer)
# Final Histogram Plot
plt.figure(figsize=(5.3, 3))  # rectangular, matches your 5.3:3 logic

sns.histplot(
    data=qwer,
    x=0,
    stat="probability",
    bins=50,
    element="step"
)

plt.axvline(x=teamscore_wt, color='r', linestyle='-')
plt.text(
    teamscore_wt - 0.03,
    0.025,
    'WT Team Strength(4.6))',
    fontsize=10,
    color='red',
    rotation=90
)

plt.xlabel('Team Strength')
plt.ylabel('Frequency')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.title('Team Strength vs Frequency', fontsize = 12)
plt.tight_layout()
plt.savefig(
    "./figures/Fig4/Fig4F_histogram.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()
