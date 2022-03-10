stk.v.12.0
WrittenBy    STK_v12.0.0

BEGIN Chain

    Name		 Chain_target_orbit2POI
    BEGIN Definition

        Object		 Place/Paris
        Object		 Satellite/Target_orbit
        Type		 Chain
        FromOperator		 Or
        FromOrder		 1
        ToOperator		 Or
        ToOrder		 1
        Recompute		 Yes
        IntervalType		 0
        ComputeIntervalStart		 0
        ComputeIntervalStop		 86400
        ComputeIntervalPtr		
        BEGIN EVENTINTERVAL
            BEGIN Interval
                Start		 16 Aug 2021 00:00:00.000000000
                Stop		 17 Aug 2021 00:00:00.000000000
            END Interval
            IntervalState		 Explicit
        END EVENTINTERVAL

        ConstConstraintsByStrands		 Yes
        UseSaveIntervalFile		 No
        SaveIntervalFile		 D:\Documents\Spacecraft-formation-control\numerical_studies\Letters_impuslive_control\Mission_design\STK\strand.int
        UseMinAngle		 No
        UseMaxAngle		 No
        UseMinLinkTime		 No
        LTDelayCriterion		 2
        TimeConvergence		 0.005
        AbsValueConvergence		 1e-14
        RelValueConvergence		 1e-08
        MaxTimeStep		 360
        MinTimeStep		 0.01
        UseLightTimeDelay		 Yes
        DetectEventsUsingSamplesOnly		 No
        UseLoadIntervalFile		 No
        BEGIN StrandObjIndexes
            StrandObj		 Place/Paris
            StrandObj		 Satellite/Target_orbit
        END StrandObjIndexes

        SaveMode		 1
        BEGIN StrandAccessesByIndex
            Strand		 0 1
            Start		  8.9830470819180791e+03
            Stop		  9.4855432094491262e+03
            Start		  1.4693047733700183e+04
            Stop		  1.5540280221261592e+04
            Start		  2.0635807692148031e+04
            Stop		  2.1533782189667083e+04
            Start		  2.6862606458199469e+04
            Stop		  2.7420023463636826e+04
            Start		  5.8511819352707840e+04
            Stop		  5.9086253422299116e+04
            Start		  6.4403296297396621e+04
            Stop		  6.5303397816996694e+04
            Start		  7.0398118629480276e+04
            Stop		  7.1240967389561149e+04
            Start		  7.6455096148341298e+04
            Stop		  7.6944750190646737e+04
        END StrandAccessesByIndex


    END Definition

    BEGIN Extensions

        BEGIN ExternData
        END ExternData

        BEGIN ADFFileData
        END ADFFileData

        BEGIN Desc
            BEGIN ShortText

            END ShortText
            BEGIN LongText

            END LongText
        END Desc

        BEGIN Crdn
        END Crdn

        BEGIN Graphics

            BEGIN Attributes

                StaticColor		 #ffffff
                AnimationColor		 #00ff00
                AnimationLineWidth		 2
                StaticLineWidth		 3

            END Attributes

            BEGIN Graphics
                ShowGfx		 On
                ShowStatic		 Off
                ShowAnimationHighlight		 On
                ShowAnimationLine		 On
                ShowLinkDirection		 Off
            END Graphics
        END Graphics

        BEGIN VO
        END VO

    END Extensions

END Chain

