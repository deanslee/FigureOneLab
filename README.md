**Figure One Lab** (**F1L**) is an initiative to help you build a computational biology (compbio) portfolio by re-enacting Figure 1 of biology papers.

**Challenge 1: Building Complex Portfolio Projects**

While many excellent university courses, online tutorials, and in-person workshops have been created to teach compbio, most of them don't help you build a portfolio that approximates the full complexity of compbio work. They often feature projects with pre-cleaned data and an expected analysis outcome, which are too simple relative to what working computational biologists grapple with daily.

F1L bridges this complexity gap. Using F1L to re-enact Figure 1 of biology papers, you can create portfolio projects that better reflect the full complexity of compbio work.

**Challenge 2: Knowing Where to Start**

So you are eager to start replicating Figure 1. But which paper do you choose? It is often hard to know where to start.

A strong compbio portfolio project consists of a broadly interesting **question/paper** (questions are usually described in papers), relevant **data**, and **methods** to analyze said data. Because the number of possible question-data-method combinations is so large, most people get decision paralysis and never start. Of those who do manage to get started, many find that they are stuck with questions that are too niche, inaccessible data, or poorly maintained methods, and then they are back to square one.

F1L shepherds you through these initial pitfalls so you don't get stuck from the start.

**Solution: F1L**

To address Challenge 1 and 2, F1L has thoughtfully curated a set of questions/papers, data, and methods to get you started:
1. **Questions/papers** are sourced from the fields of neuroscience, stem cell biology, cancer biology, and immunology. Figure 1 of these papers serve as our starting point. These papers address topics of broad appeal to both academia and industry so you are working on problems that people actually care about.
   - [Trevino_2021_Cell](https://www.sciencedirect.com/science/article/pii/S0092867421009429) Human cortex development
   - [Uzquiano_2022_Cell](https://www.sciencedirect.com/science/article/pii/S0092867422011680) Human cortical organoids
   - [Schmidt_2022_Science](https://www.science.org/doi/10.1126/science.abj4008) CRISPR screen in human T cells
   - [Kinker_2020_NatureGenetics](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8135089/) Biological programs active across commonly used cancer cell lines
   - More will be added over time.
2. **Data** are open (not behind a paywall or in a restricted database), human, and multi-omics. F1L emphasizes single-cell RNA-sequencing (scRNA-seq) data at first because it is a rich data modality ubiquitous in every field of biology and a good gateway data modality that links to many other data modalities.
3. **Methods** are well-documented Python packages. F1L does not emphasize method development because there are already plenty of powerful, popular methods that you can use to put together a brilliant portfolio. No need to re-invent the wheel.

**How It Works**

Choose a paper that piques your interest. Read it, probably several times. Skim through the supplemental figures and files as well. Use Google and ChatGPT to clarify difficult concepts. 

Then go through the Jupyter notebooks accompanying the paper. In these notebooks I demonstrate how a working computational biologist might think about reconstructing Figure 1. Run each code block and see if you can understand what was done and why. **The notebooks get you examining the data quickly but does not clarify much of the ambiquity in the data. Untangling that ambiguity is an important aspect of compbio work I have reserved for you.** Refer back to the paper frequently for the context of the data.

The notebooks only recreate a part of Figure 1. Try to reproduce the rest of Figure 1, **but don't worry if you don't produce an exact replica**. That is not the point. As you wrestle with the data behind Figure 1, you should notice nuances in the interpretation of the data that warrant further investigation. Let your curiosity and scientific instinct determine your next set of questions and analyses. This might lead you to reproduce other figures of the paper or bring in data from another paper of your own choosing. This is where your project gets exciting.

Your final product can be a little or a lot, depending on how much follow-up investigation you do.

If you need further guidance, you can email me directly at dean.sl.lee@gmail.com. I will do my best to give you pointers with the time I have.

**Intended Audience**

- Bare minimum prerequisites: 1 year of college-level biology (AP Biology also works), introductory Python, habit of googling or ChatGPT-ing to finds answers
- Anyone looking to build a strong portfolio for a compbio internship or job, in either academia or industry.
- Biologists trying to transition from wet lab to compbio work.
- Anyone with quantitative training - statistics, applied math, computer science, etc. - trying to transition to compbio work.
- Computational biologists trying to develop an expertise in new data modalities or fields of biology.
- If you are already analyzing data in an academic lab for a manuscript with the help of a mentor, you don't need F1L.
