stk.v.12.0
WrittenBy    STK_v12.0.0

BEGIN Chain

    Name		 Chain1
    BEGIN Definition

        Object		 Place/Paris
        Object		 Satellite/Satellite1
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
            StrandObj		 Satellite/Satellite1
        END StrandObjIndexes

        SaveMode		 1
        BEGIN StrandAccessesByIndex
            Strand		 0 1
            Start		  9.0193711541187749e+03
            Stop		  9.5103131351271604e+03
            Start		  1.4713107120253917e+04
            Stop		  1.5551339669014127e+04
            Start		  2.0640288613226341e+04
            Stop		  2.1529165688574125e+04
            Start		  2.6851541051829074e+04
            Stop		  2.7398307851618789e+04
            Start		  5.8421581312632814e+04
            Stop		  5.8971212543005815e+04
            Start		  6.4291621138568931e+04
            Stop		  6.5180859294716225e+04
            Start		  7.0269671537227376e+04
            Stop		  7.1107196065816010e+04
            Start		  7.6311089761647643e+04
            Stop		  7.6799903521397573e+04
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

                StaticColor		 #00ff00
                AnimationColor		 #00ffff
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

