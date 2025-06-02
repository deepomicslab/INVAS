//
// Created by caronkey on 10/5/2021.
//

#ifndef SEQMAP_WEIGHT_H
#define SEQMAP_WEIGHT_H


namespace seqGraph {
    class Weight {
    protected:
        float mCoverage;           // coverage as adjusted by LP
        float mCopyNum;             // the copy number, coverage / average coverage of one copy
        float mCopyNumOriginal;
        float mCopyNumBackup;
        float conCoverage;

        bool mIsInferred;
        float flow;
        float rev;
        float capacity;

    public:
        Weight(float aCoverage);

        ~Weight();

        float getCoverage();

        float getCopyNum();

        float getCopyNumBackup();

        void setCoverage(float aCoverage);

        void setCopyNum(float aCopyNum);

        void backup();

        void restore();

        void increaseCopyNum(float aIncrement = 1);

        void decreaseCopyNum(float aDecrement = 1);

        bool isInferred();

        void setInferred();

        void resetInferred();

        void print();

        float getFlow();
        void setFlow(float aFlow);

        float getCapacity();
        void setCapacity(float capacity);

        void setConCoverage(float conCoverage);
        float getConCoverage();
    };
}


#endif //SEQMAP_WEIGHT_H
