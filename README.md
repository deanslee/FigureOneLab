## **Figure One Lab**

### **Figure One Lab** (**F1L**) is a makerspace to help you hone computational biology skills by re-enacting Figure 1 of modern biology papers.

### **What is in F1L?**

F1L contains a curated set of questions, data, and methods that you can use to build a strong compbio project:
1. **Questions** presented in the form of papers are sourced from the fields of cancer biology, immunology, neuroscience, stem cell biology. Figure 1 of these papers serve as our starting point. These papers address topics of broad scientific appeal to both academia and industry so you are working on problems that people actually care about.
   - [Kinker_2020_NatureGenetics](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8135089/) is study of biological programs in commonly used cancer cell lines.
   - [Schmidt_2022_Science](https://www.science.org/doi/10.1126/science.abj4008) is a study of CRISPR-edited human T cells.
   - [Trevino_2021_Cell](https://www.sciencedirect.com/science/article/pii/S0092867421009429) is a study of human cortex development.
   - [Uzquiano_2022_Cell](https://www.sciencedirect.com/science/article/pii/S0092867422011680) is a study of human cortical organoids.
2. **Data** accompanying these papers are open (not behind a paywall or in a restricted database), from human tissue, and multiomic. F1L emphasizes single-cell RNA-sequencing (scRNA-seq) data at first because it is a rich data modality ubiquitous in every field of biology and a good gateway data modality that can be integrated with many other types of data.
3. **Methods** are well-documented Python packages. F1L does not emphasize method development because there are already plenty of powerful, popular methods that you can use to put together a brilliant compbio project. No need to reinvent the wheel.

### **Who Should Use F1L?**

- Bare minimum prerequisites: relentless curiosity, 4+ hours to commit a week, 1 year of college-level biology (AP Biology also works), introductory Python, habit of googling or chatting your way to answers.
- Bachelor's/master's/PhD students and postdocs in the life sciences who wish to grow/deepen their compbio skills while building a strong portfolio for a compbio internship or job.
- Scientists trying to transition from wet lab to compbio work.
- Anyone with quantitative training but not in relation to the life sciences - statistics, applied math, computer science, etc. - trying to transition to compbio work.
- If you are already analyzing data in an academic or industry setting with the help of a compbio mentor, you probably don't need F1L. Take full advantage of your mentor!

### **How to Start**

Choose one of the F1L papers that piques your interest. Read it, probably several times. Sleuth around to determine how the authors shared their data; then download it. Skim through the supplemental figures and files as well. Use Google and ChatGPT to clarify difficult concepts.

Then go through the Jupyter notebooks accompanying the paper. Note that these are not polished notebooks by any means. They are rough drafts, so you can decide how to polish them. F1L believes in doing quick-and-dirty analyses first to develop an intuition for the data before working on the final version of your code/analyses. In these notebooks I demonstrate how a working computational biologist might think about reconstructing Figure 1. Run each code block. See if you can understand what was done and why.

The notebooks only recreate a fraction of Figure 1, leaving much of the ambiguity surrounding the data for you to untangle. **Try to reproduce the rest of Figure 1, but know that the point isn't to produce an exact replica.** The point is that as you examine and wrestle with the data behind Figure 1, you should notice nuances in the interpretation of the data that warrant further investigation. Let your curiosity and scientific instinct determine your next set of questions. This might lead you to validate your ideas with other data from the same paper or with data from a new paper. This is where your project gets exciting!

Your final product might be small or large, depending on how much follow-up investigation you do.

A simple directory structure that you might use to organize your work:

<img width="325" alt="Screen Shot 2024-01-01 at 5 42 52 PM" src="https://github.com/deanslee/FigureOneLab/assets/35471368/43a129e4-fa33-469a-9206-f07d9854f071">

- FigureOneLab/trevino/data/ holds all the "raw" files, PDFs, code provided by authors, Excel spreadsheets, etc., from the paper needed to re-enact Figure 1.

- FigureOneLab/trevino/outs/ holds all the intermediate and final outputs of your analysis.

- FigureOneLab/trevino/filename1.ipynb is your first analysis.

- FigureOneLab/trevino/filename2.ipynb is your second analysis. You get the drift.

### **A Concise List of Helpful Resources**
- [3Blue1Brown](https://www.youtube.com/c/3blue1brown) for linear algebra
- [Points of Significance](https://www.nature.com/collections/qghhqm/pointsofsignificance) for statistics
- [StatQuest!!!](https://statquest.org/) for statistics
- [STAT115 by Shirley Liu](https://www.youtube.com/playlist?list=PLeB-Dlq-v6tY3QLdQBA7rwb4a7fK9mLpv) for an introduction to compbio
- [Best practices for scRNA-seq data analysis](https://www.sc-best-practices.org/preamble.html) the best resource of its kind so far
- [Glittr](https://glittr.org/?per_page=25&sort_by=stargazers&sort_direction=desc) a collection of Git repos with bioinformatics training material
