## **Figure One Lab**

### **Figure One Lab** (**F1L**) is an initiative to help you build a computational biology (compbio) portfolio by re-enacting Figure 1 of biology papers.

### **Challenge 1: Doing Sufficiently Complex Projects**

While many excellent university courses, online tutorials, and in-person workshops have been created to teach compbio, most of them don't help you build a portfolio that approximates the full complexity of compbio work. They often feature projects with pre-cleaned data and an expected analysis outcome, which are too simple relative to what working computational biologists grapple with daily.

F1L bridges this complexity gap. Using F1L to re-enact Figure 1 of biology papers, you can create portfolio projects that better reflect the full complexity of compbio work.

### **Challenge 2: Knowing Where to Start**

So you are eager to start replicating Figure 1. But which paper do you choose? It is often hard to know where to start.

A strong compbio portfolio project consists of a broadly interesting **question/paper** (questions are usually described in papers), relevant **data**, and **methods** to analyze said data. Because the number of possible question-data-method combinations is so large, most people get decision paralysis and never start. Of those who do manage to get started, many find that they are stuck with questions that are too niche, inaccessible data, or poorly maintained methods, and then they are back to square one.

F1L shepherds you through these initial pitfalls so you don't get stuck from the start.

### **Solution: F1L**

To address Challenge 1 and 2, F1L has thoughtfully curated a set of questions/papers, data, and methods to get you started:
1. **Questions/papers** are sourced from the fields of neuroscience, stem cell biology, cancer biology, and immunology. Figure 1 of these papers serve as our starting point. These papers address topics of broad appeal to both academia and industry so you are working on problems that people actually care about.
   - [Trevino_2021_Cell](https://www.sciencedirect.com/science/article/pii/S0092867421009429) Human cortex development
   - [Uzquiano_2022_Cell](https://www.sciencedirect.com/science/article/pii/S0092867422011680) Human cortical organoids
   - [Schmidt_2022_Science](https://www.science.org/doi/10.1126/science.abj4008) CRISPR screen in human T cells
   - [Kinker_2020_NatureGenetics](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8135089/) Biological programs active across commonly used cancer cell lines
   - [Hedou_2024_NatureBiotechnology](https://www.nature.com/articles/s41587-023-02033-x) Finding omic biomarkers, coming soon
   - More will be added over time.
2. **Data** accompanying these papers are open (not behind a paywall or in a restricted database), human, and multi-omics. F1L emphasizes single-cell RNA-sequencing (scRNA-seq) data at first because it is a rich data modality ubiquitous in every field of biology and a good gateway data modality that links to many other data modalities.
3. **Methods** are well-documented Python packages. F1L does not emphasize method development because there are already plenty of powerful, popular methods that you can use to put together a brilliant portfolio. No need to reinvent the wheel.

### **How It Works**

Choose one of the F1L papers that piques your interest. Read it, probably several times. Sleuth around to determine how the authors shared their data; then download it. Skim through the supplemental figures and files as well. Use Google and ChatGPT to clarify difficult concepts.

Then go through the Jupyter notebooks accompanying the paper. In these notebooks I demonstrate how a working computational biologist might think about reconstructing Figure 1. Run each code block. See if you can understand what was done and why. The notebooks only recreate a fraction of Figure 1, leaving much of the ambiguity surrounding the data for you to untangle. Try to reproduce the rest of Figure 1, but know that **the point isn't to produce an exact replica**. The point is that as you examine and wrestle with the data behind Figure 1, you should notice nuances in the interpretation of the data that warrant further investigation. Let your curiosity and scientific instinct determine your next set of questions. This might lead you to validate your ideas with other data from the same paper or with data from a new paper. This is where your project gets exciting!

Your final product might be small or large, depending on how much follow-up investigation you do.

A simple directory structure that you might use to organize your work:

<img width="325" alt="Screen Shot 2024-01-01 at 5 42 52 PM" src="https://github.com/deanslee/FigureOneLab/assets/35471368/43a129e4-fa33-469a-9206-f07d9854f071">

- FigureOneLab/trevino/data/ holds all the "raw" files, PDFs, code provided by authors, Excel spreadsheets, etc., from the paper needed to re-enact Figure 1.

- FigureOneLab/trevino/outs/ holds all the intermediate and final outputs of your analysis.

- FigureOneLab/trevino/filename1.ipynb is your first analysis.

- FigureOneLab/trevino/filename2.ipynb is your second analysis. You get the drift.

### **Intended Audience**

- Bare minimum prerequisites: 1 year of college-level biology (AP Biology also works), introductory Python, habit of googling or ChatGPT-ing to find answers
- Anyone looking to build a strong portfolio for a compbio internship or job, in either academia or industry.
- Biologists trying to transition from wet lab to compbio work.
- Anyone with quantitative training - statistics, applied math, computer science, etc. - trying to transition to compbio work.
- Computational biologists trying to develop an expertise in new data modalities or fields of biology.
- If you are already analyzing data in an academic or industry setting with the help of a compbio mentor, you probably don't need F1L. Take full advantage of your mentor!

### **A Concise List of Helpful Resources**
- [3Blue1Brown](https://www.youtube.com/c/3blue1brown) for linear algebra
- [Points of Significance](https://www.nature.com/collections/qghhqm/pointsofsignificance) for statistics
- [StatQuest!!!](https://statquest.org/) for statistics
- [STAT115 by Shirley Liu](https://www.youtube.com/playlist?list=PLeB-Dlq-v6tY3QLdQBA7rwb4a7fK9mLpv) for an introduction to compbio
- [Best practices for scRNA-seq data analysis](https://www.sc-best-practices.org/preamble.html) the best resource of its kind so far

### Mentoring

If you need further guidance with how to use F1L, you can email me directly at dean.sl.lee@gmail.com. Put "F1L" somewhere in the subject of the email. I will do my best to give you pointers with the time I have.
