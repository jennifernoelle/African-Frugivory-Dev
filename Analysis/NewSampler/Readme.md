Information for final runs

- No pseudoabsences: Pseudoabsence results did not make sense: for example heldout pseudo-absences had a higher posterior mean than overall. This likely results from a lack of good pseduo-absences in the dataset. Only around 20 observations had 3 or more repeated non-edges, and it may be that these were due to seasonality (e.g. not fruiting at the time) rather than true non-interactions. Hence, pseudo-absences will be omitted from final runs. 
- No baboons or baboon-only plants
- New sampler: There seem to be marginal benefits to the new sampler, so we will include it in the increased runs
- Number of iterations, replications etc: 
- Scenarios: Default 0/1, Defualt 0.5/1, Default 0.75/1, Exp1mod - all with and without the prior
- Additionally, do some short runs to verify that Exp1mod beats the other scenarios
- Do we want posterior networks for all scenarios? Or just the one that has the best cv results?

    
- To do: 
  - X Make sure the preliminary results can be plotted and look good
  - X If yes, then update all of the .sh files to use the same settings
  - CV RUNS: including the best expert from the preliminary results
    - X Note that old_p75 needs to be rerun because it was initially run on too few iterations
    - X Double check that all other runs had the right number of iterations
    - Status: Adding p99 and p1 scenarios for context, rerunning old models due to error in 'detect'
    - X Extract preliminary posterior networks from CV analysis and review
  - "Final runs": run the best models without CV for post network and trait matching
    - Status: 
        - newp75 first seed failed to save properly (use the "more output" results)
        - newmod1 was terminated and started over: us either it or the more output result, depending on which finishes first
        - rerunning the old models due to error in 'detect' - but using MO version
    - Save the posterior interaction matrices and heatmaps, number of new pairs
    - How different are the resulting posterior interaction matrices
    - How different are the occurrence probabilities in the new sampler vs old sampler
      - Note: in the old sampler, we sample occurrence, not occurrence probabilities - we could make it return mean of the occurrence indicator?
    - How many scenarios should Camille actually include
    - Note: new sampler code not currently returning occurrences
      - X Update new code to return running mean plant occurrence probabilities and occurrences
      - X Update old code to return running mean of plant occurrence indicators
      - X Create new batch files to run for the four best scenarios
      - Status: running. Note this does not need to be done before Camille can rerun defaunation simulations
        - Restarted old models MO version due to error in 'detect' in the original function
        - Restarted old models MO version for 1/0.75 because O_P not updating
        - Need to rerun CV now
      
  - Run even longer final runs: mixing isn't great
    - Modify code to return only running mean of L's
    - Note: mixing is better using the fixed occurrence probabilities vs expert scenarios (at least in old model)
  - Run trait matching
    - Update code

    
    
