You are an autonomous research assistant.

Your task is defined below. Execute it fully without asking for clarification.
If something is ambiguous, make reasonable assumptions and proceed (but in
written outputs, make note of the ambiguity).

--- TASK START ---

# Task 1: Replicating Paper Figures
Our work is based heavily on the work in the paper "Parameter-setting heuristic
for the quantum alternating operator ansatz", which can be found in the
`References` folder.

I would like to be able to reproduce the graphs from each figure in that paper.

When applicable, use code already present in the repo to do so. In particular,
there are already interfaces for running QAOA, running proxy QAOA, and sampling
the paper proxy in this repo.

In general, customizability should be emphasized. In other words, I would like to
be able to produce slight variants of the graphs in ways that I will describe in
more detail below.

Also emphasize simplicitly, briefness, and readability of the scripts. Use
comments liberally to make the code easy to understand, in a kind of "literate
programming" style.

For each figure, after writing the script, run it and write a markdown report
summarizing the results. Not only the mathematical results, but the results of
the programming process. E.g., what was your approach for the script, what
challenges were there, what are some future directions to go? You should
especailly alert me if it is unclear how to replicate the paper results based on
the description in the paper.

Also write a summary of the script to CLAUDE.md, such that if I make a variation
of the figures, Claude will know where to look to see existing work that can be
copied or modified.

Makie should be used for all plotting (GLMakie or CairoMakie is fine).

After making a script, run it with small values (e.g. small number of nodes in
the graph, number of samples, etc) to quickly produce graphs and make sure that
the results are acceptable (i.e. no errors, things are being drawn correctly).

I will now provide specific instructions for each figure.

Many duplicate functionalities exist in the julia and python parts of this repo.
you should prefer the Julia ones.

Also, you should timestamp each script with when you first started working on
the file and when you finished.

While working on the scripts, write your plans in research_task1_plan.md.

## Figure 1
Skip this one, as it is not a computational result.

## Figure 2
It should be easy to swap out the graph being characterized. It should also be
easy to save the graphs for all possible values of c' at once.

## Figure 3
It should be easy to swap out the graph being characterized. It should also be
able to swap out the "analytical method" used to caclulate N(c',d,c). What the
paper refers to as the "analytical method" is talking about the different
proxies available. Like with figure 2, it should also be easy to save graphs for
various values of c'. It would also be nice to be able to compare multiple
proxies in one graph. And to select which points are used to make the subplot
heatmaps.

# Figure 4
As before, it should be easy to swap out which homogeneous proxy is being used.

# Figure 5
As before, it should be easy to swap out which homogeneous proxy is being used.

# Figure 6
For this one, in addition to the parameter transfer approach and paper
homogeneous proxy, I would like to compare multiple other proxies.

# Figure 7
I would like to compare multiple proxies. Also, I think the linear ramp
scheduling is important enough that it should be added to the general API, not
just confined to this one script.

--- TASK END ---
